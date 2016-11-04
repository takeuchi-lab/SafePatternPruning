#include <cstdio>
#include <cstdlib>
#include <random>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <vector>
#include <map>
#include <ctime>
#include <string.h>

using namespace std;

struct node {
	vector<int> key;
	vector<int> x;
	double val;
};

struct feature {
	vector<int> x;
	double w;
};

struct ball {
	vector<double> center;
	double radius;
};

typedef map< vector<int>, feature > model;

int T;
int n,d;
int maxdepth;
int f;
int b;
vector< vector<int> > Z;
vector<double> Y;

void read(int argc, char **argv);
void read_data(char *filename);
node get_maxval(const vector<double> &v);
void get_maxval(const vector<double> &v, vector<int> &parentkey, const vector<int> &parentx, int depth, int j, node &maxnode);
void safeprune(const ball &ball, vector<node> &safeset);
void safeprune(const ball &ball, vector<node> &safeset, vector<int> &parentkey, const vector<int> &parentx, int depth, int j);

int main(int argc, char **argv) {
	read(argc, argv);

	double bias = 0;

	// residual
	vector<double> r(n);
	for (int i = 0; i < n; i++) r[i] = Y[i];

	node maxnode = get_maxval(r);
	double lammax = maxnode.val;
	cout << "------------------------------" << endl;
	cout << "lambda_max = " << lammax << endl;

	model model;
	model[maxnode.key] = feature{maxnode.x, 0};

	vector<double> theta(n);
	for (int i = 0; i < n; i++) theta[i] = r[i]/lammax;

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
		vector< vector<int> > key(m);
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

				double maxval = 0;
				if (iter == 1) {
					node maxnode = get_maxval(r);
					maxval = maxnode.val;
				} else {
					for (int s = 0; s < m; s++) {
						int j = index[s];
						int xnorm = x[j].size();
						double xtr = 0;
						for (int i = 0; i < xnorm; i++) {
							int idx = x[j][i];
							xtr += r[idx];
						}
						if (fabs(xtr) > maxval) maxval = fabs(xtr);
					}
				}

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
					// update model
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
				double radius = sqrt(2*gap)/lam;
				vector<node> safeset;
				int active = 0;
				if (iter == 1) { // first pruning
					map<vector<int>, feature> old;
					for (int j = 0; j < m; j++) {
						if (w[j] != 0) {
							old[key[j]] = feature{x[j], w[j]};
						}
					}

					ball B = {theta, radius};
					safeprune(B, safeset);

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

node get_maxval(const vector<double> &v) {
	node maxnode;
	maxnode.val = 0;

	vector<int> key;
	key.reserve(maxdepth);
	get_maxval(v, key, vector<int>(), 0, 0, maxnode);
	return maxnode;
}

void get_maxval(const vector<double> &v, vector<int> &parentkey, const vector<int> &parentx, int depth, int j, node &maxnode) {
	for (; j < d; j++) {
		parentkey.push_back(j);

		vector<int> x; // and
		if (depth) {
			set_intersection(Z[j].begin(), Z[j].end(), parentx.begin(), parentx.end(), back_inserter(x));
		} else {
			x = Z[j];
		}

		if (!x.size()) {
			parentkey.pop_back();
			continue;
		}

		// ここから内積計算
		double xTv = 0;
		double p = 0;
		double m = 0;
		for (int i = 0; i < (int)x.size(); i++) {
			double val = v[x[i]];
			xTv += val;
			(val > 0) ? p += val : m += val;
		}

		if (max(p,-m) < maxnode.val) {
			parentkey.pop_back();
			continue;
		}

		double score = fabs(xTv);
		if (score > maxnode.val) {
			maxnode.key = parentkey;
			maxnode.x = x;
			maxnode.val = score;
		}

		if (j < d-1 && depth+1 < maxdepth) {
			get_maxval(v, parentkey, x, depth+1, j+1, maxnode);
		}

		parentkey.pop_back();
	}
}

void safeprune(const ball &ball, vector<node> &safeset) {
	vector<int> key;
	key.reserve(maxdepth);
	safeprune(ball, safeset, key, vector<int>(), 0, 0);
}

void safeprune(const ball &ball, vector<node> &safeset, vector<int> &parentkey, const vector<int> &parentx, int depth, int j) {
	for (; j < d; j++) {
		parentkey.push_back(j);

		vector<int> x; //and
		if (depth) {
			set_intersection(Z[j].begin(), Z[j].end(), parentx.begin(), parentx.end(), back_inserter(x));
		} else {
			x = Z[j];
		}

		if (!x.size()) {
			parentkey.pop_back();
			continue;
		}

		// centerとxの内積
		double xTc = 0;
		double p = 0;
		double m = 0;
		for (int i = 0; i < (int)x.size(); i++) {
			double val = ball.center[x[i]];
			xTc += val;
			(val > 0) ? p += val : m += val;
		}

		if (max(p, -m) + ball.radius*sqrt(x.size()) < 1) {
			parentkey.pop_back();
			continue;
		}

		double score = fabs(xTc) + ball.radius*sqrt(x.size());
		if (score >= 1) {
			safeset.push_back(node{parentkey, x, score});
		}

		if (j < d-1 && depth+1 < maxdepth) {
			safeprune(ball, safeset, parentkey, x, depth+1, j+1);
		}

		parentkey.pop_back();
	}
}

void read(int argc, char **argv) {
	int i;
	T = 100;
	maxdepth = 3;
	f = 100;
	b = 0;
	char filename[1024];
	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i >= argc)
			exit(1);
		switch (argv[i-1][1]) {
		case 'T':
			T = atoi(argv[i]);
			break;
		case 'D':
			maxdepth = atoi(argv[i]);
			break;
		case 'F':
			f = atoi(argv[i]);
			break;
		case 'B':
			b = atoi(argv[i]);
			break;
		default:
			cout << "unknown option" << endl;
			exit(1);
			break;
		}
	}

	if (i >= argc) exit(1);
	strcpy(filename, argv[i]);
	read_data(filename);
}

char* readline(FILE *input) {
	int max_line_len = 1024;
	char* line = (char *)malloc(max_line_len*sizeof(char));

	int len;
	if (fgets(line,max_line_len,input) == NULL) return NULL;

	while (strrchr(line, '\n') == NULL) {
		max_line_len *= 2;
		line = (char *) realloc(line,max_line_len);
		len = (int) strlen(line);
		if (fgets(line+len, max_line_len-len, input) == NULL) break;
	}
	return line;
}

void read_data(char *filename) {
	FILE *fp = fopen(filename, "r");
	if (fp == NULL) {
		fprintf(stderr, "Cannot open input file\n");
		exit(1);
	}

	vector< vector<int> > rowZ;

	n = 0;
	d = 0;
	double Ysum = 0;
	char* line = NULL;
	while ((line = readline(fp)) != NULL) {
		rowZ.push_back(vector<int>());

		char *y = strtok(line, " \t\n");
		Y.push_back(atof(y));
		Ysum += Y[n];

		while(1) {
			char *idx = strtok(NULL, ":");
			char *val = strtok(NULL, " \t");
			if (val == NULL) break;

			int j = atoi(idx);
			rowZ[n].push_back(j);
			if (j > d) d = j;
		}
		n++;
	}
	Ysum /= n;
	for (int i = 0; i < n; i++) {
		Y[i] -= Ysum;
	}

	Z.reserve(d);
	for (int j = 0; j < d; j++) {
		Z.push_back(vector<int>());
	}

	for (int i = 0; i < n; i++) {
		for (int j = 0; j < (int)rowZ[i].size(); j++) {
			int idx = rowZ[i][j]-1;
			Z[idx].push_back(i);
		}
	}
}
