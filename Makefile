CC=g++
CFLAGS=-O3
LDFLAGS=-static

pdc: pdc.cpp PDBParser.hpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
