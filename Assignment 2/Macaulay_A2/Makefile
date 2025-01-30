# Variables
CXX = g++
CXXFLAGS = -std=c++17 -Wall
LDFLAGS =
SOURCES = gameof1d.cpp initialize.cpp time_step.cpp output.cpp
OBJECTS = $(SOURCES:.cpp=.o)
TARGET = gameof1d

# Default target
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) -o $@ $^

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $<

# Clean up generated files
clean:
	rm -f $(OBJECTS) $(TARGET)
