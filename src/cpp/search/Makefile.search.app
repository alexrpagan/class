###  BASIC PROJECT SETTINGS
APP = search
SRC = search

EWAH_ROOT=../../../ewah
TABIX_ROOT=../../../vcflib/tabixpp

LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = tabix boost_system boost_filesystem $(LIB_:%=%$(STATIC))

LIBS = $(TABIX_ROOT)/tabix.o \
       $(TABIX_ROOT)/bgzf.o \
       $(TABIX_ROOT)/knetfile.o \
       $(TABIX_ROOT)/index.o \
       $(TABIX_ROOT)/kstring.o \
       $(CMPRS_LIBS) $(DL_LIBS) $(NETWORK_LIBS) $(ORIG_LIBS) -L$(TABIX_ROOT)

CPPFLAGS = $(ORIG_CPPFLAGS) -I$(EWAH_ROOT)/headers/ -I$(TABIX_ROOT)

