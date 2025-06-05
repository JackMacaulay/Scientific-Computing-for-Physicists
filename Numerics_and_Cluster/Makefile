# Makefile for the friction module and utilities
#
# On the cluster, you will need to load at the gcc, rarray, and boost
# modules before typing 'make'
#
# To create the test model data, type "make testmodel.dat".
# To analyze the data, type "make analysis".
#
# Part of assignment 4 of PHY1610H, Winter 2025

CXX=g++
CXXFLAGS=-O3 -g -Wall -Wfatal-errors -std=c++17
LDLIBS=-lboost_program_options

all: friction.o testmodel analyze

friction.o: friction.cpp

testmodel.o: testmodel.cpp model.h
model.o: model.cpp model.h
analyze.o: analyze.cpp friction.h

testmodel: model.o testmodel.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

analyze: analyze.o friction.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(LDLIBS)

testmodel.dat: testmodel
	./testmodel -f $@ -a 0.125 -g 9.8 -z 0.0 -v 10.0 -t 16.0 -d 0.25

analysis.out: testmodel.dat analyze
	./analyze -f $< -o $@

analysis: analysis.out
	cat analysis.out

clean:
	$(RM) testmodel.o model.o friction.o analyze.o analyze testmodel analysis.out testmodel.dat

doc:
	    doxygen Doxyfile

.PHONY: clean all analysis analysis.out
