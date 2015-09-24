### prima/Makefile
###
### Copyright 2013 Elias Assmann

ifeq "$(origin FC)" "default"
	FC = gfortran
endif

F2PY=f2py --quiet --fcompiler=$(FC)

sppv.so: SpaghettiPrimavera.f90
	$(F2PY) -c -m sppv --f90flags=-fbounds-check \
		only: spaghettiprimavera sppv_data \
		      charactertorgbcolor charactertothickness \
		: SpaghettiPrimavera.f90
