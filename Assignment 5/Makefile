# Makefile for Part 1 and Part 2 of the KPP-Fisher assignment

CXX       = g++
# Note the extra -I path for cblas.h
CXXFLAGS  = -O2 -std=c++11 -I/mingw64/include/openblas

# Common libraries (math library, etc.)
LIBS      = -lm

# BLAS library (for Part 2 only).
# Also specify -L/mingw64/lib so the linker can find libopenblas.dll.a
BLASLIBS  = -L/mingw64/lib -lopenblas -lm

all: fisher_loop_explicit fisher_blas

######################################################
# Part 1 target
######################################################
fisher_loop_explicit: fisher_loop_explicit.cpp
	$(CXX) $(CXXFLAGS) fisher_loop_explicit.cpp -o fisher_loop_explicit $(LIBS)

######################################################
# Part 2 target
######################################################
fisher_blas: fisher_blas.cpp
	$(CXX) $(CXXFLAGS) fisher_blas.cpp -o fisher_blas $(BLASLIBS)

######################################################
# Cleanup
######################################################
clean:
	rm -f fisher_loop_explicit fisher_blas
