###  PATH TO A PRE-BUILT C++ TOOLKIT
include Makefile.builddir

###  LIST OF PROJECTS TO BUILD, SUITABLE FOR IN-TREE USE
include Makefile.out

###  RULE TO REGENERATE Makefile.out
%.out: %.in
	 sed -e 's,@builddir@,$$(builddir),g; s,@srcdir@,.,g' $< > $@

###  USE LOCAL TREE -- DON'T EDIT OR MOVE THESE TWO LINES !!!
MAKE_LIB = $(MAKE) -f Makefile.$${i}_lib $(MFLAGS)
MAKE_APP = $(MAKE) -f Makefile.$${i}_app $(MFLAGS)

###  DIRECT THE CHECK FRAMEWORK TO LOOK IN THE RIGHT PLACE
import_root = .
CHECK_ADD_KET += ; sed -e 's:[^ ]*\( ____ \):.\1:' check.sh.list > .csl &&     mv .csl check.sh.list

###  PUT YOUR OWN ADDITIONAL TARGETS (MAKE COMMANDS/RULES) BELOW HERE
