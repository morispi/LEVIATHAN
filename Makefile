curDir = $(shell pwd)
CC = g++
CFLAGS  = -Wall -pedantic -O3 -m64 -shared -Wl,--whole-archive -lpthread -Wl,--no-whole-archive -std=c++11 -fPIC

BAMTOOLS_INC = $(curDir)/LRez/bamtools/include/bamtools/
BAMTOOLS_LIB = $(shell readlink -f $(curDir)/LRez/bamtools/lib*)

LREZ_INC = $(curDir)/LRez/src/include/
LREZ_LIB = $(curDir)/LRez/lib/

MYVC_INC = $(curDir)/src/include/
MYVC_LIB = $(curDir)/lib/

CTPL_INC = $(curDir)/CTPL/
SSW_INC = $(curDir)/Complete-Striped-Smith-Waterman-Library/src/

LDFLAGS_BAMTOOLS = -lbamtools -L$(BAMTOOLS_LIB)
LDFLAGS_LREZ = -llrez -L$(LREZ_LIB)
LDFLAGS_MYVC = -lpthread

MAIN = src/main.o
SOURCE = src/alignmentsProcessing.o src/barcodesProcessing.o src/candidatesProcessing.o src/misc.o src/supportComputation.o src/SVValidation.o src/globalVariables.o src/genomeProcessing.o src/SVOuput.o src/help.o

EXEC = bin/LEVIATHAN

all: directories $(EXEC)

directories:
	mkdir -p bin/

$(EXEC): $(MAIN) $(SOURCE)
	$(CC) -o $(EXEC) $(MAIN) $(SOURCE) $(LDFLAGS_BAMTOOLS) -Wl,-rpath,$(BAMTOOLS_LIB) $(LDFLAGS_LREZ) -Wl,-rpath,$(LREZ_LIB) $(LDFLAGS_MYVC)

src/%.o: src/%.cpp
	$(CC) -o $@ -c $< $(CFLAGS) -I$(BAMTOOLS_INC) -I$(LREZ_INC) -I$(MYVC_INC) -I$(CTPL_INC) -I$(SSW_INC)


clean:
	rm src/*.o $(EXEC)


.PHONY: all clean directories
