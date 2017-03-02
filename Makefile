# C++ compiler and flags
CXX=g++
CXXFLAGS=-std=c++11 -Wall -Wextra -Wpedantic -Wshadow -Wnon-virtual-dtor \
	 -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual \
	 -Wconversion -Wsign-conversion -Wmisleading-indentation -Wnonnull \
	 -fopenmp

# Libraries and includes
LDLIBS=-lgsl -lgslcblas -lm -lgomp
# LDFLAGS=
# INC=

# Files
SOURCES=dbotu.cpp dbotu_opts.cpp common.cpp
OBJECTS=$(SOURCES:.cpp=.o)
EXECUTABLE=dbotu

# Debug builld
debug: CXXFLAGS += -g -DDEBUG -O0
debug: all

# Static build
debug_static: CXXFLAGS += -static -O3
debug_static: all

# Release build
release: CXXFLAGS += -DNDEBUG -O3
release: all

release_static: CXXFLAGS += -DNDEBUG -static -O3
release_static: all

# Rules
all: $(EXECUTABLE)


$(EXECUTABLE): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o $(EXECUTABLE) $^ $(LDLIBS)

.cpp.o:
	$(CXX) -c $(CXXFLAGS) $(INC) $< -o $@

.PHONY: clean
clean:
	rm $(OBJECTS) $(EXECUTABLE)
