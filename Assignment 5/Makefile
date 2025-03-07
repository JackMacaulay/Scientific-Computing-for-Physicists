CXX       = g++

CXXFLAGS  = -O2 -std=c++11 -I/mingw64/include/openblas

LIBS      = -lm

BLASLIBS  = -L/mingw64/lib -lopenblas -lm

all: fisher_loop_explicit fisher_blas fisher_sparse

fisher_loop_explicit: fisher_loop_explicit.cpp
	$(CXX) $(CXXFLAGS) fisher_loop_explicit.cpp -o fisher_loop_explicit $(LIBS)


fisher_blas: fisher_blas.cpp
	$(CXX) $(CXXFLAGS) fisher_blas.cpp -o fisher_blas $(BLASLIBS)

fisher_sparse: fisher_sparse.cpp
	$(CXX) $(CXXFLAGS) fisher_sparse.cpp -o fisher_sparse $(LIBS)


clean:
	rm -f fisher_loop_explicit fisher_blas fisher_sparse
