CXXFLAGS = -std=c++17 -O3 -Wall -g -IC:/msys64/mingw64/include
LDFLAGS = -g
TEST_LIBS = -lCatch2Main -lCatch2

all: hydrogen

hydrogen.o: hydrogen.cpp eigenval.h output.h init.h
eigenval.o: eigenval.cpp eigenval.h
output.o: output.cpp output.h
init.o: init.cpp init.h
output_c2.o: output_c2.cpp output.h
test_init.o: test_init.cpp init.h
test_eigenval.o: test_eigenval.cpp eigenval.h
test_integration.o: test_integration.cpp init.h eigenval.h output.h

hydrogen: hydrogen.o eigenval.o output.o init.o
	$(CXX) $(LDFLAGS) -o $@ $^

output_c2: output_c2.o output.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(TEST_LIBS)

test_init: test_init.o init.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(TEST_LIBS)

test_eigenval: test_eigenval.o eigenval.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(TEST_LIBS)

test_integration: test_integration.o init.o eigenval.o output.o
	$(CXX) $(LDFLAGS) -o $@ $^ $(TEST_LIBS)

run_output_c2: output_c2
	./output_c2 -s

run_test_init: test_init
	./init_test -s

run_test_eigenval: test_eigenval
	./test_eigenval -s

run_test_integration: test_integration
	./test_integration -s

test: run_output_c2 run_test_init run_test_eigenval run_test_integration

clean:
	$(RM) *.o hydrogen output_c2 test_init test_eigenval test_integration

.PHONY: all clean run_output_c2 run_test_init run_test_eigenval run_test_integration test
