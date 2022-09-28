CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=#-static

all: pdc pdd

pdc: pdc.cpp PDBParser.hpp StringTools.hpp pstream.h Superpose.hpp GeometryTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdd: pdd.cpp PDBParser.hpp StringTools.hpp pstream.h Superpose.hpp GeometryTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
