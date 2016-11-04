CXX ?= g++;
CC ?= gcc
CFLAGS ?= -Wall -std=gnu++11 -O3 -fPIC

all: train

train: ./gspan/gspan.o ./gspan/graph.o ./gspan/misc.o ./gspan/ismin.o ./gspan/dfs.o train.cc
	$(CXX) $(CFLAGS) -o train train.cc ./gspan/gspan.o ./gspan/graph.o ./gspan/misc.o ./gspan/ismin.o ./gspan/dfs.o
	
./gspan/gspan.o: ./gspan/gspan.cc ./gspan/gspan.h ./gspan/tree.hh
	$(CXX) $(CFLAGS) -c -o ./gspan/gspan.o ./gspan/gspan.cc
	
./gspan/graph.o: ./gspan/graph.cc ./gspan/gspan.h ./gspan/tree.hh
	$(CXX) $(CFLAGS) -c -o ./gspan/graph.o ./gspan/graph.cc
	
./gspan/misc.o: ./gspan/misc.cc ./gspan/gspan.h
	$(CXX) $(CFLAGS) -c -o ./gspan/misc.o ./gspan/misc.cc
	
./gspan/ismin.o: ./gspan/ismin.cc ./gspan/gspan.h
	$(CXX) $(CFLAGS) -c -o ./gspan/ismin.o ./gspan/ismin.cc
	
./gspan/dfs.o: ./gspan/dfs.cc ./gspan/gspan.h
	$(CXX) $(CFLAGS) -c -o ./gspan/dfs.o ./gspan/dfs.cc
	
clean:
	rm -f *~ ./gspan/gspan.o ./gspan/graph.o ./gspan/misc.o ./gspan/ismin.o ./gspan/dfs.o train