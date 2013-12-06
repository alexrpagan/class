###  BASIC PROJECT SETTINGS
APP = search
SRC = search timer

EWAH_ROOT=../../../ewah
VCFLIB_ROOT=../../../vcflib
TABIX_ROOT=$(VCFLIB_ROOT)/tabixpp
FASTAHACK_ROOT=$(VCFLIB_ROOT)/fastahack

LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = tabix boost_system boost_filesystem $(LIB_:%=%$(STATIC))

LIBS = $(TABIX_ROOT)/tabix.o \
       $(TABIX_ROOT)/bgzf.o \
       $(TABIX_ROOT)/knetfile.o \
       $(TABIX_ROOT)/index.o \
       $(TABIX_ROOT)/kstring.o \
       $(FASTAHACK_ROOT)/Fasta.o \
       $(FASTAHACK_ROOT)/split.o \
       $(CMPRS_LIBS) $(DL_LIBS) $(NETWORK_LIBS) $(ORIG_LIBS) -L$(TABIX_ROOT)

CPPFLAGS = $(ORIG_CPPFLAGS) -I$(EWAH_ROOT)/headers/ -I$(TABIX_ROOT) -I$(FASTAHACK_ROOT)

