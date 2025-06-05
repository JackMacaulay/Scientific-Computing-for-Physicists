# Compiler settings
CXX      = g++
CXXFLAGS = -O2 -std=c++17 -Wall -Wextra

# Source file and target executable name
SRC      = macaulay_a7.cpp
TARGET   = macaulay_a7

# Default target: build the executable.
all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

# Run the simulation with default parameters.
# Parameters: M=200, N=200, p=0.7, g=2.0, K=25, S=20
run: $(TARGET)
	./$(TARGET) 200 200 0.7 2.0 25 20

# Clean up generated files.
clean:
	rm -f $(TARGET)
