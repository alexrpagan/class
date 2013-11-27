# Compiler
CC=g++
C=gcc

# Compiler flags
CFLAGS=-O3 -D_FILE_OFFSET_BITS=64
#CFLAGS=-O3 -static -D VERBOSE_DEBUG  # enables verbose debugging via --debug2

SOURCE_ROOT=src/cpp
VCFLIB_ROOT=vcflib

LIBS = -L./ -L$(VCFLIB_ROOT)/tabixpp/ -ltabix -lz -lm
INCLUDE = -I$(VCFLIB_ROOT)/src -I$(VCFLIB_ROOT)/

all: bin/genbitsets

static:
	$(MAKE) CFLAGS="$(CFLAGS) -static" all

debug:
	$(MAKE) CFLAGS="$(CFLAGS) -D VERBOSE_DEBUG -g -rdynamic" all

profiling:
	$(MAKE) CFLAGS="$(CFLAGS) -g" all

gprof:
	$(MAKE) CFLAGS="$(CFLAGS) -pg" all

.PHONY: all

OBJECTS=split.o \
		$(VCFLIB_ROOT)/tabixpp/tabix.o \
		$(VCFLIB_ROOT)/tabixpp/bgzf.o \
		$(VCFLIB_ROOT)/smithwaterman/SmithWatermanGotoh.o \
		$(VCFLIB_ROOT)/smithwaterman/disorder.c \
		$(VCFLIB_ROOT)/smithwaterman/LeftAlign.o \
		$(VCFLIB_ROOT)/smithwaterman/Repeats.o \
		$(VCFLIB_ROOT)/smithwaterman/IndelAllele.o \
		Variant.o \

# executables

genbitsets bin/genbitsets: genbitsets.o $(OBJECTS) $(HEADERS)
	$(CC) $(CFLAGS) $(INCLUDE) genbitsets.o $(OBJECTS) -o bin/genbitsets $(LIBS)

# objects

genbitsets.o: $(SOURCE_ROOT)/genbitsets.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SOURCE_ROOT)/genbitsets.cpp

clean:
	rm -f *.o
	rm -f bin/*
