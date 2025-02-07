# Makefile
CXXFLAGS=-std=c++17 -O3 -Wall -g
LDFLAGS=-g
all: hydrogen
hydrogen.o: hydrogen.cpp eigenval.h output.h init.h
eigenval.o: eigenval.cpp eigenval.h
output.o: output.cpp output.h
init.o: init.cpp init.h
output_c2.o: output_c2.cpp output.h
hydrogen: hydrogen.o eigenval.o output.o init.o
	$(CXX) $(LDFLAGS) -o $@ $^
output_c2: output_c2.o output.o
	$(CXX) $(LDFLAGS) -o $@ $^ -lCatch2Main -lCatch2
run_output_c2: ./output_c2
	./output_c2 -s
clean:
	$(RM) hydrogen.o eigenval.o output.o init.o output_c2.o
.PHONY: run_output_c2 all clean
