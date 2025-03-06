CXX = g++
CXXFLAGS = -O2 -std=c++11
LIBS = -lm

all: fisher_loop_explicit

fisher_loop_explicit: fisher_loop_explicit.cpp
	$(CXX) $(CXXFLAGS) fisher_loop_explicit.cpp -o fisher_loop_explicit $(LIBS)

clean:
	rm -f fisher_loop_explicit
