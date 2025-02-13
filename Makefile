CXX = g++
CXXFLAGS = -O2 -Wall -fPIC -std=c++23 -Wno-deprecated-declarations -D_GLIBCXX_USE_CXX11_ABI=0

# ROOT flags and libs
ROOTCFLAGS    = $(shell root-config --cflags)
ROOTLIBS      = $(shell root-config --libs)
ROOTGLIBS     = $(shell root-config --glibs)
ROOT_LIBDIR   = $(shell root-config --libdir)

# Pythia8 configuration
PYTHIA8_INCLUDE = $(PYTHIA8)/include
PYTHIA8_LIB = $(PYTHIA8)/lib
PYTHIA8_FLAGS = -I$(PYTHIA8_INCLUDE)
PYTHIA8_LIBS = -L$(PYTHIA8_LIB) -Wl,-rpath,$(PYTHIA8_LIB) -lpythia8

# FastJet configuration
FASTJET_INCLUDE = $(FASTJET)/include
FASTJET_LIB = $(FASTJET)/lib
FASTJET_FLAGS = -I$(FASTJET_INCLUDE)
FASTJET_LIBS = -L$(FASTJET_LIB) -lfastjet

# Additional flags for FastJet compatibility
EXTRA_FLAGS = -D_BACKWARD_BACKWARD_WARNING_H -I.

# Combine all flags
ALL_FLAGS = $(CXXFLAGS) $(ROOTCFLAGS) $(PYTHIA8_FLAGS) $(FASTJET_FLAGS) $(EXTRA_FLAGS)
LIBS = $(ROOTLIBS) $(PYTHIA8_LIBS) $(FASTJET_LIBS)

# Source files and targets
LIB_SRCS = TaggingUtilities.cpp
LIB_OBJS = $(LIB_SRCS:.cpp=.o)
LIB_TARGET = libTaggingUtilities.so

# bjet_analysis targets
ANALYSIS_TARGET = bjet_analysis
ANALYSIS_SRCS = bjet_analysis.C

# Add rpath for dynamic libraries (removing duplicates)
UNIQUE_RPATHS = $(sort $(PYTHIA8_LIB) $(FASTJET_LIB) $(ROOT_LIBDIR) .)
RPATH = $(addprefix -Wl,-rpath,$(UNIQUE_RPATHS))

all: $(LIB_TARGET) $(ANALYSIS_TARGET)

$(LIB_TARGET): $(LIB_OBJS)
	$(CXX) -shared -o $@ $(LIB_OBJS) $(LIBS) $(RPATH)

$(ANALYSIS_TARGET): $(ANALYSIS_SRCS) $(LIB_TARGET)
	$(CXX) $(ALL_FLAGS) -o $@ $(ANALYSIS_SRCS) -L. -lTaggingUtilities $(LIBS) $(RPATH) -DSTANDALONE_MODE

%.o: %.cpp
	$(CXX) $(ALL_FLAGS) -c $< -o $@

clean:
	rm -f $(LIB_TARGET) $(ANALYSIS_TARGET)
	rm -f *.o
	rm -f *~
	rm -f *.d

.PHONY: all clean