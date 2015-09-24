### prima/Makefile
###
### Copyright 2013 Elias Assmann

ifeq "$(origin FC)" "default"
	FC = gnu95
endif

ifneq "$(findstring gnu95,$(FC))" ""
	FFLAGS := -O3 $(FFLAGS)
else ifneq "$(findstring intel,$(FC))" ""
	FFLAGS := -O2 $(FFLAGS)
endif


F2PY ?= f2py --quiet --fcompiler=$(FC)

sppv.so: SpaghettiPrimavera.f90
	$(F2PY) -c -m sppv --f90flags="$(FFLAGS)" \
		only: spaghettiprimavera sppv_data \
		      charactertorgbcolor charactertothickness \
		: SpaghettiPrimavera.f90

clean:
	rm -f sppv.so *.mod
