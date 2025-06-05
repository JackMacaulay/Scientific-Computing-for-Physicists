CXX       = g++
CXXFLAGS  = -O2 -std=c++11 -I/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/openblascore/0.3.24/include
LIBS      = -lm
BLASLIBS  = -L/cvmfs/soft.computecanada.ca/easybuild/software/2023/x86-64-v4/Compiler/gcccore/openblascore/0.3.24/lib -lopenblas -lm

all: fisher_loop_explicit fisher_blas

fisher_loop_explicit: fisher_loop_explicit.cpp
	$(CXX) $(CXXFLAGS) fisher_loop_explicit.cpp -o fisher_loop_explicit $(LIBS)

fisher_blas: fisher_blas.cpp
	$(CXX) $(CXXFLAGS) fisher_blas.cpp -o fisher_blas $(BLASLIBS)

clean:
	rm -f fisher_loop_explicit fisher_blas
