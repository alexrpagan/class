###  BASIC PROJECT SETTINGS
APP = search
SRC = search

EWAH_ROOT=../../../ewah

LIB_ = $(BLAST_INPUT_LIBS) $(BLAST_LIBS) $(OBJMGR_LIBS)
LIB = boost_system boost_filesystem $(LIB_:%=%$(STATIC))

LIBS = $(CMPRS_LIBS) $(DL_LIBS) $(NETWORK_LIBS) $(ORIG_LIBS)

CPPFLAGS = $(ORIG_CPPFLAGS) -I$(EWAH_ROOT)/headers/

