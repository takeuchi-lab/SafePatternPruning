CXX ?= g++;
CC ?= gcc
CFLAGS ?= -Wall -std=gnu++11 -O3 -fPIC

all: train

train: train.cc
	$(CXX) $(CFLAGS) -o train train.cc
	
clean:
	rm -f *~ train