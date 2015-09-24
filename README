
                             prima.py 0.2
			     ************

This is a wrapper script for Maurits Haverkort's `SpaghettiPrimavera.f90´.
It produces the same kind of plots but should provide a simple
interface.  Like `SpaghettiPrimavera.f90´, it needs a `.qtl_band´ and
a `.klist_band´ file; furthermore it reads a `.struct´ if present.

See `prima.py -h´ and `prima.py´ (executed in a suitable `case´
directory) for usage hints.

By  way of  example,  suppose we  are  interested in  the  t2g and  eg
contributions of the first inequivalent atom in the struct.  A typical
command line might be

   prima.py 1:D-t2g:5,0,0:3 1:D-eg:0,0,5:3 -o prima.ps

if  that atom  has  `ISPLIT=2´.   If instead  it  has `ISPLIT=8´,  the
command would read:

   prima.py 1:DXY+DXZ+DYZ:5,0,0:3 1:DX2Y2+DZ2:0,0,5:3 -o prima.ps

which is an abbreviated form of

   prima.py 1:DXY:5,0,0:3 1:DXZ:5,0,0:3 1:DYZ:5,0,0:3 \
            1:DX2Y2:0,0,5:3 1:DZ2:0,0,5:3 -o prima.ps


Installation
============

prima.py needs Python 2, NumPy, F2PY (which is now part of NumPy), and
a modern Fortran compiler.  I tested it with the following versions:

Python 	 2.7.3
NumPy  	 1.6.1
F2PY   	 2
gfortran 4.6.3


Comments and Suggestions
========================

about prima.py should go to <elias.assmann@gmail.com>.

The author of the original SpaghettiPrimavera.f90 is
Maurits W. Haverkort <Maurits@physics.ubc.ca>.
