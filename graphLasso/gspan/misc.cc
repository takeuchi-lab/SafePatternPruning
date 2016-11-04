#include "gspan.h"

namespace GSPAN {

const RMPath& DFSCode::buildRMPath()
{
	rmpath.clear(); // 現在のDFSCodeのrmpath

	int old_from = -1;

	for (int i = size() -1; i >= 0; --i) {
		if ((*this)[i].from < (*this)[i].to && // 辺がforward
				(rmpath.empty() || old_from == (*this)[i].to)) { // 1つ前の辺の始点が,今の辺の終点(つながっている)
			rmpath.push_back(i);
			old_from = (*this)[i].from;
		}
	}

	return rmpath;
}

void History::build (Graph &graph, PDFS *e)
{
	// first build history
	clear ();
	edge.clear ();   edge.resize   (graph.edge_size());
	vertex.clear (); vertex.resize (graph.size());

	if (e) { // eがnullでない
		push_back (e->edge);
		edge  [e->edge->id] = vertex[e->edge->from] =  vertex[e->edge->to]   = 1; // 辺番号 = 辺の始点 = 辺の終点 = 1
		for (PDFS *p = e->prev; p ; p = p->prev) { // pがnullでない間
			push_back (p->edge);
			edge  [p->edge->id] = vertex[p->edge->from] =  vertex[p->edge->to]   = 1;
		}
		std::reverse (begin(), end());
	}
}

bool get_forward_rmpath (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear ();
	int tolabel = graph[e->to].label; // eの終点のlabel

	for (Vertex::edge_iterator it = graph[e->from].edge.begin(); it != graph[e->from].edge.end(); ++it) { // グラフのeの始点から出てる辺について
		int tolabel2 = graph[it->to].label; // その辺の終点ラベル
		if (e->to == it->to || minlabel > tolabel2 || history.hasVertex (it->to)) continue; // eの終点 == その辺の終点 or minlabel > tolabel2 or historyがその終点を持っている なら飛ばす
		if (e->elabel < it->elabel || (e->elabel == it->elabel && tolabel <= tolabel2)) // eのlabel < その辺のlabel or (eのlabel == その辺のlabel and eの終点label <= その辺の終点label)
			result.push_back (&(*it));
	}

	return (! result.empty());
}

bool get_forward_pure (Graph &graph, Edge *e, int minlabel, History& history, EdgeList &result)
{
	result.clear ();

	for (Vertex::edge_iterator it = graph[e->to].edge.begin(); it != graph[e->to].edge.end(); ++it) { // グラフのeの終点から出てる辺について
		if (minlabel > graph[it->to].label || history.hasVertex (it->to)) continue; // minlabel > その辺の終点ラベル or historyがその終点を持っている なら飛ばす
		result.push_back (&(*it));
	}

	return (! result.empty());
}

bool get_forward_root(Graph &g, Vertex &v, EdgeList &result) {
	result.clear();
	for (Vertex::edge_iterator it = v.edge.begin(); it != v.edge.end(); ++it) // vから出る辺すべて
		if (v.label <= g[it->to].label) result.push_back(&(*it)); //vのラベルより大きいラベルに進むとき

	return (!result.empty()); // 空でないなら1
}

Edge *get_backward(Graph &graph, Edge* e1, Edge *e2, History& history)
{
	if (e1 == e2) return 0;

	for (Vertex::edge_iterator it = graph[e2->to].edge.begin(); it != graph[e2->to].edge.end(); ++it) {
		if (history.hasEdge(it->id)) continue;
	      if ( (it->to == e1->from) &&
		   ( (e1->elabel < it->elabel)  ||
		     (e1->elabel == it->elabel) && (graph[e1->to].label <= graph[e2->to].label)
			 ) ) return &(*it);

	}

	return 0;
}

}
