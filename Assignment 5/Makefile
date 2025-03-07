# Makefile for Part 1, Part 2, and Part 3 of the KPP-Fisher assignment

CXX       = g++
# Use the extra -I flag for the BLAS header in Part 2.
CXXFLAGS  = -O2 -std=c++11 -I/mingw64/include/openblas

# Common libraries (e.g., math library)
LIBS      = -lm

# BLAS library (for Part 2 only).
# Specify -L/mingw64/lib so the linker can find libopenblas.dll.a.
BLASLIBS  = -L/mingw64/lib -lopenblas -lm

all: fisher_loop_explicit fisher_blas fisher_sparse

######################################################
# Part 1 target: Full matrix with explicit triple loop.
######################################################
fisher_loop_explicit: fisher_loop_explicit.cpp
	$(CXX) $(CXXFLAGS) fisher_loop_explicit.cpp -o fisher_loop_explicit $(LIBS)

######################################################
# Part 2 target: Full matrix with BLAS routine.
######################################################
fisher_blas: fisher_blas.cpp
	$(CXX) $(CXXFLAGS) fisher_blas.cpp -o fisher_blas $(BLASLIBS)

######################################################
# Part 3 target: Sparse matrix approach.
######################################################
fisher_sparse: fisher_sparse.cpp
	$(CXX) $(CXXFLAGS) fisher_sparse.cpp -o fisher_sparse $(LIBS)

######################################################
# Cleanup
######################################################
clean:
	rm -f fisher_loop_explicit fisher_blas fisher_sparse
