#include "gspan.h"
#include <cstring>
#include <string>
#include <iterator>
#include <set>

namespace GSPAN {

// 現在のDFSCodeのグラフgを生成
bool DFSCode::toGraph(Graph &g)
{
	g.clear();

	for (DFSCode::iterator it = begin(); it != end(); ++it) {
		g.resize(std::max(it->from, it->to) + 1);

		if (it->fromlabel != -1) g[it->from].label = it->fromlabel;
		if (it->tolabel != -1) g[it->to].label = it->tolabel;
		g[it->from].push(it->from, it->to, it->elabel);
		g[it->to].push(it->to, it->from, it->elabel);
	}

	g.buildEdge();

	return true;
}

std::ostream &DFSCode::write (std::ostream &os)
{
	if (size() == 0) return os;

	os << "(" << (*this)[0].fromlabel << ") " << (*this)[0].elabel << " (0f" << (*this)[0].tolabel << ")";

	for (unsigned int i = 1; i < size(); ++i) {
		if ((*this)[i].from < (*this)[i].to) {
			os << " " << (*this)[i].elabel << " (" << (*this)[i].from << "f" << (*this)[i].tolabel << ")";
		} else {
			os << " " << (*this)[i].elabel << " (b" << (*this)[i].to << ")";
		}
	}

	return os;
}

}
