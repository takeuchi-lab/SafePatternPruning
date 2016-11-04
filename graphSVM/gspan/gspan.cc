#include "gspan.h"
#include <strstream>
#include <algorithm>
#include <functional>
#include <numeric>
#include <fstream>
#include <cstdlib>
#include <math.h>
#include <ctime>

using namespace std;

#define LEARN 0;
#define CLASSIFY 1

namespace GSPAN {

void gSpan::read(int argc, char **argv) {
	int i;
	for (i = 1; i < argc; i++) {
		if (argv[i][0] != '-') break;
		if (++i >= argc)
			exit(1);
		switch (argv[i-1][1]) {
		case 'T':
			T = atoi(argv[i]);
			break;
		case 'F':
			f = atoi(argv[i]);
			break;
		case 'S':
			minsup = atoi(argv[i]);
			break;
		case 'D':
			maxpat = atoi(argv[i]);
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

	ifstream fp(argv[i]);
	if (!fp) {
		fprintf(stderr, "Cannot open input file\n");
		exit(1);
	}
	read_data(fp);
}

std::istream &gSpan::read_data(std::istream &is)
{
	Graph g;
	while (true) {
		g.read (is);
		if (g.empty()) break;
		transaction.push_back (g);
	}
	return is;
}

void gSpan::get_maxval(void) {
	rule.gain = 0.0;

	std::vector<tree<TNODE>::pre_order_iterator> mount;
	mount.clear();

	tree<TNODE>::pre_order_iterator pre=tr.begin();
	while(pre != tr.end()) {

		DFS_CODE.clear();
		for (std::vector<DFS>::iterator dit = pre->dfscode.begin(); dit != pre->dfscode.end(); dit++) {
			DFS_CODE.push_back(*dit);
		}

		bool p = false;
		p = calc_score(pre->projected);

		if ((pre->id == 1) && !p) {
			mount.push_back(pre);
		}

		pre++;
	}

	for (std::vector<tree<TNODE>::pre_order_iterator>::iterator it = mount.begin(); it != mount.end(); ++it) {

		tree<TNODE>::pre_order_iterator pre2;
		pre2 = *it;
		DFS_CODE.clear();
		for (std::vector<DFS>::iterator dit = pre2->dfscode.begin(); dit != pre2->dfscode.end(); dit++) {
			DFS_CODE.push_back(*dit);
		}
		get_maxval(pre2->projected, pre2);
	}
}

void gSpan::get_maxval(Projected &projected, tree<TNODE>::iterator &tnode)
{
	if (!is_min()) {
		return;
	}

	if (maxpat == (int)DFS_CODE.size()) {
		return;
	}

	int support = 0;
	unsigned int oid = 0xffffffff;
	for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
		if (oid != it->id)
			++support;
	}
	if (support < minsup) {
		return;
	}

	bool p = false;
	p = calc_score(projected);

	if (p && (tnode->id == 1)) {
		return;
	}

	tree<TNODE>::iterator child;

	if (!p && (tnode->id == 1)) {
		tnode->id = 0;
		child = tnode;
	} else {
		TNODE tn;

		tn.id = 0;
		tn.dfscode.clear();
		for (std::vector<DFS>::iterator it = DFS_CODE.begin(); it != DFS_CODE.end(); ++it) {
			tn.dfscode.push_back(*it);
		}

		tn.projected.clear();
		for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
			tn.projected.push(it->id, it->edge, it->prev);
		}

		if (tnode != tr.begin())
			child = tr.append_child(tnode,tn);
		else
			child = tr.insert(tnode,tn);

		if (p) {
			child->id = 1;
			return;
		}
	}

	const RMPath &rmpath = DFS_CODE.buildRMPath();
	int minlabel         = DFS_CODE[0].fromlabel;
	int maxtoc           = DFS_CODE[rmpath[0]].to;

	Projected_map3 new_fwd_root;
	Projected_map2 new_bck_root;
	EdgeList edges;

	for (unsigned int n = 0; n < child->projected.size(); ++n) {

		unsigned int id = child->projected[n].id;

		PDFS *cur = &(child->projected[n]);

		History history (transaction[id], cur);

		// backward
		for (int i = (int)rmpath.size()-1; i >= 1; --i) {
			Edge *e = get_backward (transaction[id], history[rmpath[i]], history[rmpath[0]], history);
			if (e) new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
		}

		// pure forward
		if (get_forward_pure (transaction[id], history[rmpath[0]], minlabel, history, edges)) {
			for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it) {
				new_fwd_root[maxtoc][(*it)->elabel][transaction[id][(*it)->to].label].push (id, *it, cur);
			}
		}

		// backtracked forward
		for (int i = 0; i < (int)rmpath.size(); ++i) {
			if (get_forward_rmpath (transaction[id], history[rmpath[i]], minlabel, history, edges)) {
				for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it) {
					new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][transaction[id][(*it)->to].label].push (id, *it, cur);
				}
			}
		}
	}

	// backward (recursive)
	for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
		for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
			DFS_CODE.push (maxtoc, to->first, -1, elabel->first, -1);
			get_maxval(elabel->second, child);
			DFS_CODE.pop();
		}
	}

	// forward
	for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
		for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				DFS_CODE.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
				get_maxval(tolabel->second, child);
				DFS_CODE.pop ();
			}
		}
	}
	return;
}

bool gSpan::calc_score(Projected &projected) {
	double xTr = 0;
	double p = 0;
	double m = 0;

	int support = 0;
	unsigned int oid = 0xffffffff;
	for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
		if (oid != it->id) {
			++support;
			if (r[it->id] > 0) {
				double val = Y[it->id]*r[it->id];
				xTr += val;
				(val > 0) ? p += val : m += val;
			}
		}
		oid = it->id;
	}

	if (support < minsup || max(p, -m) <= rule.gain) return true;

	double score = fabs(xTr);
	if (score > rule.gain) {
		rule.gain = score;
		rule.size = DFS_CODE.size();
		rule.dfscode = DFS_CODE;

		rule.loc.clear ();
		unsigned int oid = 0xffffffff;
		for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
			if (oid != it->id) {
				rule.loc.push_back(it->id);
			}
			oid = it->id;
		}
	}
	return false;
}


void gSpan::safeprune(vector<node> &safeset) {
	std::vector<tree<TNODE>::pre_order_iterator> mount;
	mount.clear();

	tree<TNODE>::pre_order_iterator pre=tr.begin();
	while(pre != tr.end()) {

		DFS_CODE.clear();
		for (std::vector<DFS>::iterator dit = pre->dfscode.begin(); dit != pre->dfscode.end(); dit++) {
			DFS_CODE.push_back(*dit);
		}

		bool p = false;
		p = calc_score_spp(pre->projected, safeset);

		if ((pre->id == 1) && !p) {
			mount.push_back(pre);
		}

		pre++;
	}

	for (std::vector<tree<TNODE>::pre_order_iterator>::iterator it = mount.begin(); it != mount.end(); ++it) {

		tree<TNODE>::pre_order_iterator pre2;
		pre2 = *it;
		DFS_CODE.clear();
		for (std::vector<DFS>::iterator dit = pre2->dfscode.begin(); dit != pre2->dfscode.end(); dit++) {
			DFS_CODE.push_back(*dit);
		}
		safeprune(pre2->projected, pre2, safeset);
	}
}

void gSpan::safeprune(Projected &projected, tree<TNODE>::iterator &tnode, vector<node> &safeset)
{
	if (!is_min()) {
		return;
	}

	if (maxpat == (int)DFS_CODE.size()) {
		return;
	}

	int support = 0;
	unsigned int oid = 0xffffffff;
	for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
		if (oid != it->id)
			++support;
	}
	if (support < minsup) {
		return;
	}

	bool p = false;
	p = calc_score_spp(projected, safeset);

	if (p && (tnode->id == 1)) {
		return;
	}

	tree<TNODE>::iterator child;

	if (!p && (tnode->id == 1)) {
		tnode->id = 0;
		child = tnode;
	} else {
		TNODE tn;

		tn.id = 0;
		tn.dfscode.clear();
		for (std::vector<DFS>::iterator it = DFS_CODE.begin(); it != DFS_CODE.end(); ++it) {
			tn.dfscode.push_back(*it);
		}

		tn.projected.clear();
		for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
			tn.projected.push(it->id, it->edge, it->prev);
		}

		if (tnode != tr.begin())
			child = tr.append_child(tnode,tn);
		else
			child = tr.insert(tnode,tn);

		if (p) {
			child->id = 1;
			return;
		}
	}

	const RMPath &rmpath = DFS_CODE.buildRMPath();
	int minlabel         = DFS_CODE[0].fromlabel;
	int maxtoc           = DFS_CODE[rmpath[0]].to;

	Projected_map3 new_fwd_root;
	Projected_map2 new_bck_root;
	EdgeList edges;

	for (unsigned int n = 0; n < child->projected.size(); ++n) {

		unsigned int id = child->projected[n].id;

		PDFS *cur = &(child->projected[n]);

		History history (transaction[id], cur);

		// backward
		for (int i = (int)rmpath.size()-1; i >= 1; --i) {
			Edge *e = get_backward (transaction[id], history[rmpath[i]], history[rmpath[0]], history);
			if (e) new_bck_root[DFS_CODE[rmpath[i]].from][e->elabel].push (id, e, cur);
		}

		// pure forward
		if (get_forward_pure (transaction[id], history[rmpath[0]], minlabel, history, edges)) {
			for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it) {
				new_fwd_root[maxtoc][(*it)->elabel][transaction[id][(*it)->to].label].push (id, *it, cur);
			}
		}

		// backtracked forward
		for (int i = 0; i < (int)rmpath.size(); ++i) {
			if (get_forward_rmpath (transaction[id], history[rmpath[i]], minlabel, history, edges)) {
				for (EdgeList::iterator it = edges.begin(); it != edges.end();  ++it) {
					new_fwd_root[DFS_CODE[rmpath[i]].from][(*it)->elabel][transaction[id][(*it)->to].label].push (id, *it, cur);
				}
			}
		}
	}

	// backward (recursive)
	for (Projected_iterator2 to = new_bck_root.begin(); to != new_bck_root.end(); ++to) {
		for (Projected_iterator1 elabel = to->second.begin(); elabel != to->second.end(); ++elabel) {
			DFS_CODE.push (maxtoc, to->first, -1, elabel->first, -1);
			safeprune(elabel->second, child, safeset);
			DFS_CODE.pop();
		}
	}

	// forward
	for (Projected_riterator3 from = new_fwd_root.rbegin(); from != new_fwd_root.rend(); ++from) {
		for (Projected_iterator2 elabel = from->second.begin(); elabel != from->second.end(); ++elabel) {
			for (Projected_iterator1 tolabel = elabel->second.begin();
					tolabel != elabel->second.end(); ++tolabel) {
				DFS_CODE.push (from->first, maxtoc+1, -1, elabel->first, tolabel->first);
				safeprune(tolabel->second, child, safeset);
				DFS_CODE.pop ();
			}
		}
	}
	return;
}

bool gSpan::calc_score_spp(Projected &projected, vector<node> &safeset) {
	double xTc = 0;
	double p = 0;
	double m = 0;

	int support = 0;
	unsigned int oid = 0xffffffff;
	for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
		if (oid != it->id) {
			++support;
			double val = Y[it->id]*theta[it->id];
			xTc += val;
			(val > 0) ? p += val : m += val;
		}
		oid = it->id;
	}

	if (support < minsup || max(p,-m)+radius*sqrt(support) < 1) return true;

	double score = fabs(xTc) + radius*sqrt(support);
	if (score >= 1) {
		vector<int> x;
		x.reserve(support);
		unsigned int oid = 0xffffffff;
		for (Projected::iterator it = projected.begin(); it != projected.end(); ++it) {
			if (oid != it->id) {
				x.push_back(it->id);
			}
			oid = it->id;
		}
		safeset.push_back(node{DFS_CODE, x, score});
	}
	return false;
}

} // namespace GSPAN
