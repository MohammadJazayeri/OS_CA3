# Compiler and flags
CXX := g++
CXXFLAGS := -std=c++11 -Wall -pthread

# Output executable
TARGET := VoiceFilters.out

# Source files and object files
SRCS := main_parallel.cpp filters_p.cpp read_write.cpp
OBJS := $(SRCS:.cpp=.o)

# Default rule to build the program
all: $(TARGET)

# Link object files to create the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# Compile source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Clean up generated files
clean:
	rm -f $(TARGET) $(OBJS)

# Rebuild the project
rebuild: clean all

.PHONY: all clean rebuild
