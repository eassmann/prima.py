### prima/Makefile
###
### Copyright 2013-2015 Elias Assmann

svn-rev := r$(lastword '$version:  $')

VERSION := $(svn-rev)

SHELL := /bin/bash

ifeq "$(origin FC)" "default"
	FC = gnu95
endif

ifneq "$(findstring gnu95,$(FC))" ""
	FFLAGS := -O3 -Wall $(FFLAGS)
else ifneq "$(findstring intel,$(FC))" ""
	FFLAGS := -O2 $(FFLAGS)
endif


F2PY ?= f2py --quiet --fcompiler=$(FC)

all: sppv.so prima.py

.PHONY: prima.py Makefile
prima.py Makefile:
	perl -i -pe 's/\$$version:\$$/\$$version: 1.0\$$/' $@

sppv.so: SpaghettiPrimavera.f90
	$(F2PY) -c -m sppv --f90flags="$(FFLAGS)" \
		only: spaghettiprimavera sppv_data \
		      charactertorgbcolor charactertothickness \
		: SpaghettiPrimavera.f90

clean:
	rm -f sppv.so *.mod *__genmod.f90 *.pyc

dist: dir = prima-$(VERSION)
dist: clean
	mkdir $(dir); \
	cd $(dir); \
	ln -s -t. ../{Makefile,NEWS,prima.py,prima_sample,README,README.spaghettiprimavera,SpaghettiPrimavera.f90,TODO,webcolors.py}

	tar -zchf $(dir).tar.gz $(dir)
	rm -rf $(dir)
