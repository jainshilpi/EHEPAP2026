# Makefile

# Compiler
CXX := g++ -arch x86_64

# Compiler flags
CXXFLAGS := -std=c++11 -Wall -Wextra -I/Users/shilpi2015/geant4/geant4-v11.2.0/source/externals/clhep/include/

# Linker flags
LDFLAGS := -L/Users/shilpi2015/geant4/geant4-v11.2.0-install/lib/ -lG4clhep


# ROOT flags and libraries
ROOTCFLAGS := $(shell root-config --cflags)
ROOTLIBS := $(shell root-config --libs)

$(info ROOTCFLAGS = $(ROOTCFLAGS))
$(info ROOTLIBS = $(ROOTLIBS))

# Target executable
TARGET := serc19_ecal_clustering

# Source files
SRCS := serc19_ecal_clustering.C

# Object files
OBJS := $(SRCS:.C=.o)

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -o $@ $^ $(LDFLAGS) $(ROOTLIBS)

%.o: %.C
	$(CXX) $(CXXFLAGS) $(ROOTCFLAGS) -c -o $@ $<

clean:
	rm -f $(TARGET) $(OBJS)
