DISTRIB=$(PWD)
include $(DISTRIB)/config.mk
OPTIMIZE_FLAGS=-O3 -DNDEBUG -fopenmp
DEBUG_FLAGS=-g -fopenmp

install: 
	mkdir -p bin
	mkdir -p lib

# OPTIMIZED VERSION
	mkdir -p $(PLATFORM)-$(CC)
	+make $@ -C $(PLATFORM)-$(CC)  -f $(DISTRIB)/src/Makefile DISTRIB=$(DISTRIB) DEBUG="" FLAGS="$(OPTIMIZE_FLAGS)"

# DEBUG VERSION
	mkdir -p $(PLATFORM)-$(CC)-debug
	+make $@ -C $(PLATFORM)-$(CC)-debug -f $(DISTRIB)/src/Makefile DISTRIB=$(DISTRIB) DEBUG=-debug FLAGS="$(DEBUG_FLAGS)"


# OPTIMIZED VERSION
	mkdir -p $(PLATFORM)-$(CC)
	+make $@ -C $(PLATFORM)-$(CC)  -f $(DISTRIB)/app/Makefile DISTRIB=$(DISTRIB) DEBUG="" FLAGS="$(OPTIMIZE_FLAGS)"

# DEBUG VERSION
	mkdir -p $(PLATFORM)-$(CC)-debug
	+make $@ -C $(PLATFORM)-$(CC)-debug -f $(DISTRIB)/app/Makefile DISTRIB=$(DISTRIB) DEBUG=-debug FLAGS="$(DEBUG_FLAGS)"


help:
	@echo ""
	@echo ""
	@echo "make install"
	@echo "     build optimized and debug programs"
	@echo ""
	@echo "make build_doc"
	@echo "     build code documentation using Doxygen"
	@echo "     the  main page of the generated documentation is "
	@echo "     doc/html/html/index.html"
	@echo ""
	@echo "make clean"
	@echo "     clean object files"
	@echo ""
	@echo "make realclean"
	@echo "     clean object files and Doxygen documentation"
	@echo ""
	@echo ""

build_doc:
	cd doc;doxygen Doxyfile

clean_doc:
	rm -rf doc/html

cleandistrib:
	rm -rf $(PLATFORM)-$(CC)
	rm -rf $(PLATFORM)-$(CC)-debug
	rm -f  lib/*.a
	rm -f  lib/*~ *~ src/*~
	rm -f  bin/*

clean:cleandistrib

realclean:cleandistrib clean_doc


