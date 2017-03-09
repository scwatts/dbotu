# Compiler flags
AM_CXXFLAGS=-std=c++11 -Wall -Wextra -fopenmp
AM_FFLAGS=-Wall -Wextra

# Compile static binary
if STATIC_COMPILE
AM_CXXFLAGS += -static
endif


# Declaring files
bin_PROGRAMS = dbotu
dbotu_SOURCES = dbotu.cpp dbotu_opts.cpp common.cpp
