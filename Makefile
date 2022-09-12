CC=g++
CFLAGS=-O3
LDFLAGS=-static

all: pdc pdd

pdc: pdc.cpp PDBParser.hpp StringTools.hpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

pdd: pdd.cpp PDBParser.hpp StringTools.hpp pstream.h
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
