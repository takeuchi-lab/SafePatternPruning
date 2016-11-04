#include "gspan.h"
#include <cstring>
#include <string>
#include <iterator>
#include <strstream>
#include <set>

namespace GSPAN {


template <class T, class Iterator>
void tokenize (const char *str, Iterator iterator)
{
	std::istrstream is (str, std::strlen(str));
	std::copy (std::istream_iterator <T> (is), std::istream_iterator <T> (), iterator);
}

// edgeにidを割り振る
void Graph::buildEdge ()
{
	char buf[512];
	std::map <std::string, unsigned int> tmp;

	unsigned int id = 0;
	for (int from = 0; from < (int)size (); ++from) { // 全頂点について, その点から出る辺をit
		for (Vertex::edge_iterator it = (*this)[from].edge.begin ();
				it != (*this)[from].edge.end (); ++it) {
			if (from <= it->to) std::sprintf (buf, "%d %d %d", from,     it->to,   it->elabel);
			else                std::sprintf (buf, "%d %d %d", it->to,   from,     it->elabel);
			if (tmp.find (buf) == tmp.end()) { // 見つからなかったら
				it->id = id;
				tmp[buf] = id;
				++id;
			} else {
				it->id = tmp[buf];
			}
		}
	}

	edge_size_ = id;
}

std::istream &Graph::read (std::istream &is)
{
	std::vector <std::string> result;
	char line[1024];

	clear ();

	while (true) {

		unsigned int pos = is.tellg ();
		if (! is.getline (line, 1024)) {
			break;
		}

		result.clear ();
		tokenize<std::string>(line, std::back_inserter (result));

		if (result.empty()) {
			// do nothing
		} else if (result[0] == "t") {
			if (! empty()) { // use as delimiter
				is.seekg (pos, std::ios_base::beg);
				break;
			} else {
				y = atof(result[3].c_str());
			}
		} else if (result[0] == "v" && result.size() >= 3) {
			unsigned int id    = atoi (result[1].c_str());
			this->resize (id + 1);
			(*this)[id].label = atoi (result[2].c_str());
		} else if (result[0] == "e" && result.size() >= 4) {
			int from   = atoi (result[1].c_str());
			int to     = atoi (result[2].c_str());
			int elabel = atoi (result[3].c_str());

			if ((int)size () <= from || (int)size () <= to) {
				std::cerr << "Fromat Error:  define vertex lists before edges" << std::endl;
				exit (-1);
			}

			(*this)[from].push (from, to,   elabel); // 始点のedgeに登録
			(*this)[to].push   (to,   from, elabel); // 終点のedgeに登録
		}
	}

	buildEdge ();

	return is;
}



}
