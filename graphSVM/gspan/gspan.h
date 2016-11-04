#include <iostream>
#include <map>
#include <vector>
#include <set>
#include <algorithm>
#include <string>
#include <cmath>
#include <list>
#include <set>
#include <cstdio>
#include <string.h>
#include "tree.hh"

using namespace std;

namespace GSPAN {

struct Edge {
	int from;
	int to;
	int elabel;
	unsigned int id;
	Edge(): from(0), to(0), elabel(0), id(0) {};
};

class Vertex
{
public:
	typedef std::vector<Edge>::iterator edge_iterator;

	int label;
	std::vector<Edge> edge;

	void push (int from, int to, int elabel)
	{
		edge.resize(edge.size()+1);
		edge[edge.size()-1].from = from;
		edge[edge.size()-1].to = to;
		edge[edge.size()-1].elabel = elabel;
		return;
	}
};

class Graph: public std::vector<Vertex> {
private:
	unsigned int edge_size_;
public:
	double y;
	unsigned int edge_size() {return edge_size_; }
	unsigned int vertex_size() {return (unsigned int)size(); }
	void buildEdge();
	std::istream &read(std::istream &);
	std::ostream &write(std::ostream &);
};

class DFS {
public:
	int from;
	int to;
	int fromlabel;
	int elabel;
	int tolabel;
	friend bool operator != (const DFS &d1, const DFS &d2) {return (! (d1 == d2));}
	friend bool operator == (const DFS &d1, const DFS &d2) {
		return (d1.from == d2.from
				&& d1.to == d2.to
				&& d1.fromlabel == d2.fromlabel
				&& d1.elabel == d2.elabel
				&& d1.tolabel == d2.tolabel);
	}
	friend bool operator < (const DFS &d1, const DFS &d2) {
		std::vector<int> v1{d1.from,d1.to,d1.fromlabel,d1.elabel,d1.tolabel};
		std::vector<int> v2{d2.from,d2.to,d2.fromlabel,d2.elabel,d2.tolabel};
		return (v1 < v2);

	}
	DFS(): from(0), to(0), fromlabel(0), elabel(0), tolabel(0) {};
};

typedef std::vector<int> RMPath;

struct DFSCode: public std::vector<DFS> {
private:
	RMPath rmpath;
public:
	const RMPath& buildRMPath();
	bool toGraph(Graph &);
	void push(int from, int to, int fromlabel, int elabel, int tolabel) {
		resize(size() + 1);
		DFS &d = (*this)[size()-1];
		d.from = from; d.to = to; d.fromlabel = fromlabel; d.elabel = elabel; d.tolabel = tolabel;
	}
	void pop() { resize(size()-1); }
	std::ostream &write(std::ostream &);
};

class Projected;

struct PDFS {
	unsigned int id;
	Edge *edge;
	PDFS *prev;
	PDFS(): id(0), edge(0), prev(0) {};
};

class History: public std::vector<Edge*>
{
private:
	std::vector<int> edge;
	std::vector<int> vertex;

public:
	bool hasEdge(unsigned int id) {return (bool)edge[id];}
	bool hasVertex(unsigned int id) {return (bool)vertex[id];}
	void build (Graph &, PDFS *);
	History() {};
	History(Graph& g, PDFS *p) { build(g, p); }
};

// misc_function
typedef std::vector<Edge*> EdgeList;

bool get_forward_pure(Graph &, Edge *, int, History &, EdgeList &);
bool get_forward_rmpath(Graph &, Edge *, int, History &, EdgeList &);
bool get_forward_root(Graph &, Vertex &, EdgeList &);
Edge *get_backward(Graph &, Edge *, Edge *, History &);

class Rule {
public:
	std::string dfs;
	double gain;
	unsigned int size;
	std::vector<int> loc;
	DFSCode dfscode;
	friend bool operator < (const Rule &r1, const Rule &r2) {
		return fabs(r1.gain) < fabs(r2.gain);
	}
};

class Projected: public std::vector<PDFS> {
public:
	void push(int id, Edge *edge, PDFS *prev) {
		resize(size()+1);
		PDFS &d = (*this)[size()-1];
		d.id = id; d.edge = edge; d.prev = prev;
	}
};

struct TNODE {
	unsigned int id;
	DFSCode dfscode;
	Projected projected;
};

class gSpan {

private:

	typedef std::map<int, std::map<int, std::map<int, Projected> > > Projected_map3;
	typedef std::map<int, std::map<int, Projected> > Projected_map2;
	typedef std::map<int, Projected> Projected_map1;
	typedef std::map<int, std::map<int, std::map<int, Projected> > >::iterator Projected_iterator3;
	typedef std::map<int, std::map<int, std::map<int, Projected> > >::reverse_iterator Projected_riterator3;
	typedef std::map<int, std::map<int, Projected> >::iterator Projected_iterator2;
	typedef std::map<int, Projected>::iterator Projected_iterator1;

	std::vector<Graph> transaction;
	DFSCode DFS_CODE;
	DFSCode DFS_CODE_IS_MIN;
	Graph GRAPH_IS_MIN;

	Rule rule;
	tree<TNODE> tr;

	bool is_min();
	bool project_is_min(Projected &);

public:

	struct node {
		DFSCode key;
		vector<int> x;
		double val;
	};

	struct feature {
		vector<int> x;
		double w;
	};

	map<DFSCode,feature> model;

	int T;
	int n;
	int f;
	int b;
	int minsup;
	int maxpat;
	vector<double> Y;
	vector<double> r;
	vector<double> theta;

	double lammax;
	double radius;

	gSpan(): T(100), f(100), b(0), minsup(1), maxpat(5){};

	int train(void);
	void read(int argc, char **argv);
	istream &read_data(std::istream &);

	void get_maxval(void);
	void get_maxval(Projected &projected, tree<TNODE>::iterator &tnode);
	bool calc_score(Projected &projected);

	void safeprune(vector<node> &safeset);
	void safeprune(Projected &projected, tree<TNODE>::iterator &tnode, vector<node> &safeset);
	bool calc_score_spp(Projected &projected, vector<node> &safeset);
};


};
