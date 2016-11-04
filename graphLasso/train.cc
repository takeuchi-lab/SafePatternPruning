#include "./gspan/gspan.h"
#include <cstdio>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <ctime>

using namespace std;

int main(int argc, char **argv) {

	GSPAN::gSpan gspan;
	gspan.read(argc, argv);
	gspan.train();

	return 0;
}


namespace GSPAN{

int gSpan::train(void) {
	n = transaction.size();

	// response
	Y.reserve(n);
	double Ysum = 0;
	for (int i = 0; i < n; i++) {
		Y.push_back(transaction[i].y);
		Ysum += Y[i];
	}
	Ysum /= n;
	for (int i = 0; i < n; i++) {
		Y[i] -= Ysum;
	}

	double bias = 0;

	// residual
	r.reserve(n);
	for (int i = 0; i < n; i++) r.push_back(Y[i]);

	// Initial Tree
	tree<TNODE>::iterator top = tr.begin();

	Projected_map3 root;
	EdgeList edges;
	DFS_CODE.clear();

	for (unsigned int id = 0; id < transaction.size(); ++id) {
		Graph &g = transaction[id];
		for (unsigned int from = 0; from < g.size(); ++from) {
			if (get_forward_root(g, g[from], edges)) {
				for (EdgeList::iterator it = edges.begin(); it != edges.end(); ++it)
					root[g[from].label][(*it)->elabel][g[(*it)->to].label].push(id,*it,0);
			}
		}
	}

	// INITIAL SETTING of the tree
	if (tr.size() == 0) {
		DFS_CODE.clear();
		for (Projected_iterator3 fromlabel = root.begin(); fromlabel != root.end(); ++fromlabel) {
			for (Projected_iterator2 elabel = fromlabel->second.begin(); elabel != fromlabel->second.end(); ++elabel) {
				for (Projected_iterator1 tolabel = elabel->second.begin(); tolabel != elabel->second.end(); ++tolabel) {
					DFS_CODE.push(0, 1, fromlabel->first, elabel->first, tolabel->first);
					TNODE tn;
					tn.id = 1;
					for (std::vector<DFS>::iterator it = DFS_CODE.begin(); it != DFS_CODE.end(); it++) {
						tn.dfscode.push_back(*it);
					}

					tn.projected.clear();
					for (Projected::iterator it = tolabel->second.begin(); it != tolabel->second.end(); ++it) {
						tn.projected.push(it->id, it->edge, it->prev);
					}
					tr.insert(top,tn);
					DFS_CODE.pop();
				}
			}
		}
	}

	get_maxval();
	lammax = rule.gain;
	cout << "------------------------------" << endl;
	cout << "lambda_max = " << lammax << endl;

	model[rule.dfscode] = feature{rule.loc, 0};

	theta.reserve(n);
	for (int i = 0; i < n; i++) theta.push_back(r[i]/lammax);

	vector<double> lambda(T-1);
	for (int t = 0; t < T-1; t++) {
		lambda[t] = lammax*pow(10,-2.0*(t+1)/(T-1));
	}

	//
	// START: solution path algorithm
	//
	for (int t = 0; t < T-1; t++) {
		double lam = lambda[t];
		cout << "------------------------------" << endl;
		printf("[%d] lambda: %.9f (%.9f)\n", t, lam, log10(lam/lammax));

		// warm start
		int m = model.size();
		vector<int> index(m);
		vector<double> w(m);
		vector< vector<int> > x(m);
		vector<DFSCode> key(m);
		int j = 0;
		for (auto it = model.begin(); it != model.end(); it++) {
			index[j] = j;
			w[j] = it->second.w;
			x[j] = it->second.x;
			key[j] = it->first;
			j++;
		}

		// START: pre_solve
		double P_old = 10000000;
		for (int iter = 0; iter <= 1000000; iter++) {

			for (int j = 0; j < m; j++) {
				int i = j + rand()%(m-j);
				swap(index[i], index[j]);
			}

			double L1norm = 0;
			for (int s = 0; s < m; s++) {
				int j = index[s];

				double xTr = 0;
				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] += w[j];
					xTr += r[idx];
				}

				if (xTr > lam) {
					w[j] = (xTr-lam)/x[j].size();
				} else if (xTr < -lam) {
					w[j] = (xTr+lam)/x[j].size();
				} else {
					w[j] = 0;
				}

				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] -= w[j];
				}
				L1norm += fabs(w[j]);
			}

			// bias
			if (b) {
				double tmp = 0;
				for (int i = 0; i < n; i++) {
					r[i] += bias;
					tmp += r[i];
				}
				bias = tmp/n;
				for (int i = 0; i < n; i++) {
					r[i] -= bias;
				}
			}

			double loss = 0;
			for (int i = 0; i < n; i++) {
				loss += r[i]*r[i];
			}
			double P_new = 0.5*loss + lam*L1norm;
			if (fabs((P_old-P_new)/P_old) < 1e-8) break;
			P_old = P_new;
		}
		// END: pre_solve

		double D = 0;
		for (int i = 0; i < n; i++) {
			D += 0.5*Y[i]*Y[i] - 0.5*lam*lam*(theta[i]-Y[i]/lam)*(theta[i]-Y[i]/lam);
		}

		//
		// START: coodinate descent
		//
		for (int iter = 1; iter <= 1000000; iter++) {
			// START: dual
			if (iter % f == 1) {
				double loss = 0;
				double yTr = 0;
				for (int i = 0; i < n; i++) {
					loss += r[i]*r[i];
					yTr += Y[i]*r[i];
				}

				// L1norm
				double L1norm = 0;
				for (int s = 0; s < m; s++) {
					int j = index[s];
					L1norm += fabs(w[j]);
				}

				get_maxval();
				double maxval = rule.gain;

				// dual feasible solution
				double alpha = min(max(yTr/(lam*loss), -1/maxval), 1/maxval);

				double P = 0.5*loss + lam*+L1norm;
				double D_new = -0.5*lam*lam*alpha*alpha*loss + lam*alpha*yTr;
				if (D < D_new) {
					D = D_new;
					for (int i = 0; i < n; i++) {
						theta[i] = alpha*r[i];
					}
				}

				double gap = P-D;
				if (gap/P < 1e-8) {
					model.clear();
					int active = 0;
					for (int s = 0; s < m; s++) {
						int j = index[s];
						if (w[j] != 0) {
							active++;
							model[key[j]] = feature{x[j], w[j]};
						}
					}
					printf("[iter %4d] primal: %.9f, dual: %.9f, gap(relative): %.9f, active: %d\n", iter, P, D, (P-D)/P, active);
					break;
				}

				// START: safe pruning
				radius = sqrt(2*gap)/lam;
				vector<node> safeset;
				int active = 0;
				if (iter == 1) { // first pruning
					map<DFSCode, feature> old;
					for (int j = 0; j < m; j++) {
						if (w[j] != 0) {
							old[key[j]] = feature{x[j], w[j]};
						}
					}

					safeprune(safeset);

					m = (int)safeset.size();
					index.resize(m);
					w.resize(m);
					x.resize(m);
					key.resize(m);
					for (int j = 0; j < m; j++) {
						auto flag = old.find(safeset[j].key);
						if (flag != old.end()) {
							w[j] = flag->second.w;
						} else {
							w[j] = 0;
						}
						if (w[j] != 0) active++;
						x[j] = safeset[j].x;
						key[j] = safeset[j].key;
						index[j] = j;
					}
				} else { // dynamic screening
					for (int s = 0; s < m; s++) {
						int j = index[s];
						int xnorm = x[j].size();
						double xtc = 0;
						for (int i = 0; i < xnorm; i++) {
							int idx = x[j][i];
							xtc += theta[idx];
						}
						if (w[j] != 0) active++;
						if (fabs(xtc)+radius*sqrt(xnorm) < 1) {
							w[j] = 0;
							m--;
							swap(index[s], index[m]);
							s--;
						}
					}
				}
				printf("[iter %4d] primal: %.9f, dual: %.9f, gap(relative): %.9f, safeset: %d, active: %d\n", iter, P, D, (P-D)/P, m, active);
				// END: safe pruning
			}
			// END: dual

			for (int j = 0; j < m; j++) {
				int i = j + rand()%(m-j);
				swap(index[i], index[j]);
			}

			// START: update wj
			for (int s = 0; s < m; s++) {
				int j = index[s];

				double xTr = 0;
				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] += w[j];
					xTr += r[idx];
				}

				if (xTr > lam) {
					w[j] = (xTr-lam)/x[j].size();
				} else if (xTr < -lam) {
					w[j] = (xTr+lam)/x[j].size();
				} else {
					w[j] = 0;
				}

				for (int k = 0; k < (int)x[j].size(); k++) {
					int idx = x[j][k];
					r[idx] -= w[j];
				}
			}
			// END: update wj

			// bias
			if (b) {
				double tmp = 0;
				for (int i = 0; i < n; i++) {
					r[i] += bias;
					tmp += r[i];
				}
				bias = tmp/n;
				for (int i = 0; i < n; i++) {
					r[i] -= bias;
				}
			}
		}
		//
		// END: coodinate descent
		//
	}
	//
	// END: solution path algorithm
	//
	return 0;
}

}
