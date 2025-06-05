# Makefile for macaulay_a6

CXX = g++
CXXFLAGS = -std=c++11 -O2
LDFLAGS = -lfftw3 -lm
TARGET = macaulay_a6
SRC = macaulay_a6.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC) $(LDFLAGS)

clean:
	rm -f $(TARGET)
