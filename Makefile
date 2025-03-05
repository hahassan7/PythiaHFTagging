# Compiler and flags
CXX = g++
CXXFLAGS = -Wall -Wextra -O2 -std=c++20 -g

# ROOT configuration
ROOT_CFLAGS = $(shell root-config --cflags)
ROOT_LIBS = $(shell root-config --libs)
ROOT_GLIBS = $(shell root-config --glibs)

# Pythia8 configuration
PYTHIA8_DIR = $(PYTHIA8)
PYTHIA8_INCLUDE = -I$(PYTHIA8_DIR)/include
PYTHIA8_LIBS = -L$(PYTHIA8_DIR)/lib -lpythia8

# FastJet configuration - use exact output from fastjet-config
FASTJET_INCLUDE = $(shell fastjet-config --cxxflags)
FASTJET_LIBS = $(shell fastjet-config --libs)

# Combined flags
ALL_INCLUDES = -I. $(ROOT_CFLAGS) $(PYTHIA8_INCLUDE) $(FASTJET_INCLUDE)
# Put FastJet first in the linking order
ALL_LIBS = $(FASTJET_LIBS) $(ROOT_GLIBS) $(PYTHIA8_LIBS) -pthread -rdynamic

# Target executable
TARGET = bjet_analysis

# Object files
OBJECTS = bjet_analysis.o TaggingUtilities.o

# Default target
all: $(TARGET)

# Main target
$(TARGET): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $(ALL_INCLUDES) $(OBJECTS) -o $@ $(ALL_LIBS) -Wl,--no-as-needed

# Explicit compilation rules for each source file
bjet_analysis.o: bjet_analysis.cpp
	$(CXX) $(CXXFLAGS) $(ALL_INCLUDES) -c $< -o $@

TaggingUtilities.o: TaggingUtilities.cpp
	$(CXX) $(CXXFLAGS) $(ALL_INCLUDES) -c $< -o $@

# Clean (excluding source files)
clean:
	rm -f $(TARGET) *.o *~ *.d
	@echo "Cleaning object files and executable (preserving source files)"

# Very clean (including generated data)
distclean: clean
	rm -f *.log
	rm -f *.pdf
	rm -f *.eps

# Run the analysis
run: $(TARGET)
	./$(TARGET)

# Debug build
debug: CXXFLAGS += -g -DDEBUG
debug: clean $(TARGET)

# Help target
help:
	@echo "Available targets for Pythia Lambda_c analysis:"
	@echo "  all              - Build the executable"
	@echo "  clean            - Remove object files and executable"
	@echo "  distclean        - Remove all generated files"
	@echo "  run              - Run the analysis"
	@echo "  debug            - Build with debug symbols"
	@echo "  help             - Show this help message"

# Dependencies
TaggingUtilities.o: TaggingUtilities.cpp
bjet_analysis.o: bjet_analysis.cpp

.PHONY: all clean distclean run debug help
