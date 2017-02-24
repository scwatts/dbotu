# C++ compiler and flags
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	 -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual \
	 -Wconversion -Wsign-conversion -Wmisleading-indentation -Wnonnull \
	 -g -O0

# Libraries and includes
# LDLIBS=
# LDFLAGS=
# INC=

# Files
SOURCES=dbotu.cpp dbotu_opts.cpp common.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=dbotu

# Rules
all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $(LDLIBS) -o $(EXECUTABLE) $^

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $< -o $@

.PHONY: clean
clean:
	rm $(OBJECTS) $(EXECUTABLE)
