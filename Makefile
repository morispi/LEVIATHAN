PREFIX ?= /usr/local
BINDIR = $(PREFIX)/bin
LIBDIR = $(PREFIX)/lib

BUILD_PREFIX = $(shell pwd)
BUILD_BINDIR = $(BUILD_PREFIX)/bin
BUILD_LIBDIR = $(BUILD_PREFIX)/lib
LREZ_LIBDIR = $(BUILD_PREFIX)/LRez/lib/

CXX ?= g++
CXXFLAGS += -Wall -pedantic -O3 -m64 -std=c++11 -fPIC
LDFLAGS += -L$(LREZ_LIBDIR) -Wl,-rpath,$(LREZ_LIBDIR) -Wl,-rpath,$(BUILD_LIBDIR)

BAMTOOLS_LIB_PREFIX = lrez_
BAMTOOLS_LIB = $(BUILD_PREFIX)/LRez/lib$(BAMTOOLS_LIB_PREFIX)bamtools$(SHLIB_EXT)

BAMTOOLS_INC = $(BUILD_PREFIX)/LRez/include/bamtools/
LREZ_INC = ./LRez/src/include/
LEVIATHAN_INC = ./src/include/
CTPL_INC = ./CTPL/
SSW_INC = ./Complete-Striped-Smith-Waterman-Library/src/

LIBS_LREZ = -llrez
LIBS_BAMTOOLS = -l$(BAMTOOLS_LIB_PREFIX)bamtools
LIBS_LEVIATHAN = -lpthread

MAIN = src/main.o
SOURCE = src/alignmentsProcessing.o src/barcodesProcessing.o src/candidatesProcessing.o src/misc.o src/supportComputation.o src/SVValidation.o src/globalVariables.o src/genomeProcessing.o src/SVOutput.o src/help.o

EXEC = $(BUILD_BINDIR)/LEVIATHAN

all: $(EXEC)

directories:
	mkdir -p $(BUILD_BINDIR)

$(EXEC): $(MAIN) $(SOURCE) directories
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $(EXEC) $(MAIN) $(SOURCE) $(LIBS_BAMTOOLS) $(LIBS_LREZ) $(LIBS_LEVIATHAN)

src/%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) $(CPPFLAGS) -I$(BAMTOOLS_INC) -I$(LREZ_INC) -I$(LEVIATHAN_INC) -I$(CTPL_INC) -I$(SSW_INC) -o $@ -c $<

install:
	mkdir -p $(DESTDIR)$(BINDIR)
	cp -a $(EXEC) $(DESTDIR)$(BINDIR)/

clean:
	rm src/*.o $(EXEC)


.PHONY: all clean directories
