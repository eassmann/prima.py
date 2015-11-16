
#                            prima.py 0.3


USAGE
-----

This is a wrapper script for [Maurits Haverkort's `SpaghettiPrimavera.f90`]
(http://www.cpfs.mpg.de/haverkort/spaghetti_primavera).
It produces the same kind of plots but should provide a simpler
interface.  Like `SpaghettiPrimavera.f90`, it needs a`$case.klist_band`
and a corresponding `$case.qtl` file; it also reads `$case.struct` if
present.

If it exists, `$case.qtl[up|dn]_band` is preferred over `$case.qtl[up|dn]`.
(Wien2k does not use `.qtl_band` but it is often convenient to rename
a `.qtl` file intended for band structures.)

See `prima.py -h` and `prima.py` (executed in a suitable [Wien2k]
(http://wien2k.at) directory) for usage hints.


### An Example

Suppose we are interested in the t2g and eg contributions of the first
inequivalent atom in the struct.  A typical command line might be
```
   prima.py 1:D-t2g:5,0,0:3 1:D-eg:0,0,5:3 -o prima.ps
```
if  that atom  has  `ISPLIT=2`.   If instead  it  has `ISPLIT=8`,  the
command would read:
```
   prima.py 1:DXY++DXZ++DYZ:5,0,0:3 1:DX2Y2++DZ2:0,0,5:3 -o prima.ps`
```
which is an abbreviated form of
```
   prima.py 1:DXY:5,0,0:3 1:DXZ:5,0,0:3 1:DYZ:5,0,0:3 \
            1:DX2Y2:0,0,5:3 1:DZ2:0,0,5:3 -o prima.ps
```
Each of the non-option arguments names an atom, an orbital, a color, a
thickness increment, and a description for the (optional) legend;
these tokens are separated by colons `:`, and each of them may be
empty.


### Specification of Atoms and Orbitals

The ampersand `&` combines two orbitals or two sites.  Since it must
be escaped in most shells, a double plus sign `++` may be used instead
(a single `+` appears in some orbital descriptions in `.qtl` files).

For both the atom and the orbital specification, you can use a number
(in the sequence of `$case.struct` or `$case.qtl`) or a name (labels
as given in the same files).  Furthermore, the following special atom
names are accepted:

 * `all`  → all atoms
 * `istl` → the interstitial
 * `Z$N`  → all atoms with atomic number `$N`
 * `$Sy`  → the atomic symbol expands to all atoms of that element

The atomic symbol applies only if you have several inequivalent atoms
of the same element.  As long as you follow the Wien2k convention for
the struct labels (atomic symbol in the first two places, a number in
the third), there can be no ambiguity with these special cases.

All atom and orbital specifications are case insensitive.


### Colors and Normalization

Colors can be named in three ways: as a comma-separated triple `r,g,b`
of red, green, blue values (nominally ∈ [0, 1]), in HTML-type
hexadecimal notation `#RRGGBB`, or one of the [147 color names defined
in SVG](http://www.w3.org/TR/SVG11/types.html#ColorKeywords).

The QTL values from Wien2k, when the interstitial is included, are
normalized to 1, but how much of this weight is contained in a
particular muffin-tin sphere depends on a number of things.  This
leads to a “normalization problem”; since there is no general
solution, `prima.py` simply multiplies the colors and thickness
increments by the QTL values and adds them.  Therefore, it can be
useful to have r,g,b > 1.  To make this easier, the color may be
prefixed by a factor: `$f × $color`.

In general, colored spaghetti are not a quantitative tool, and care
must be taken in their interpretation, especially concerning the
normalization and mixing of colors.


### Mixing spins

The options `-m|--mix-spins` and `-j|--join-spins` may be used to
combine data from two ‘qtl’ files (normally `up` and `dn` spin) in one
plot.  While `-m` adds the character from the two files together and
imposes the energies to be the same, `-j` draws separate lines.

The spin can be selected by appending `↑` / `↓` / `↑↓` (or `-up` /
`-dn` / `-updn`) to a regular atom specification, e.g. `1↑`, `all-dn`.

The file names default to `$case.qtl{up,dn}_band` and may be set by
giving `--qtl` twice.


### Configuration Files

The option `-C|--config-file` may be used to read options from a file
in addition to the command line.  Each line in the config file(s)
should contain one option and, if applicable, its argument, separated
by whitespace.  Any option valid on the command line may be used,
leading dashes may be omitted.  For instance, if, by default, you want
to suppress the minor ticks, add a legend, and output to
`$case_prima.ps`, then put the following into `~/.prima`
```
   minor-ticks	0
   write-legend
   out-suffix   _prima.ps
```
and add
```
   alias prima.py="prima.py -C ~/.prima"
```
to your `~/.bashrc` or equivalent.  Beside a global config file, it
can also be useful to have per-directory configuration.  This could be
implemented by
```
  alias prima.py='prima.py $(test -e .prima && echo "-C .prima")'
```
Parsing of the config files is primitive.  Lines are skipped if they
are empty or start with a hash `#`.  Leading whitespace and whitespace
between option and argument is stripped, but trailing whitespace is
preserved.


INSTALLATION
------------

To download `prima.py`, go to the [GitHub page](https://github.com/eassmann/prima.py)
and download from `master` for the cutting-edge version, or get the
latest [release](https://github.com/eassmann/prima.py/releases) if you
prefer more stable code.

`prima.py` needs Python 2, NumPy, F2PY (which is now part of NumPy), and
a modern Fortran compiler.  I tested it with the following versions:

 * [Python 2](https://python.org/)                2.7.9
 * [NumPy](http://numpy.org/)                     1.8.2
 * [F2PY](https://sysbio.ioc.ee/projects/f2py2e/) 2
 * [gfortran](https://gcc.gnu.org/fortran/)       4.9.2

Running `make` in the `prima` directory should produce the file
`sppv.so`.  Python must be able to find this file, and the file
`webcolors.py`.  This is easiest if you create a symlink to prima.py
in some directory in your $PATH, or add the `prima` directory to your
$PATH.


MORE INFORMATION
----------------

The [wiki](https://github.com/eassmann/prima.py/wiki) has a set of
[example plots](https://github.com/eassmann/prima.py/wiki/Examples)
made with prima.py.  Questions regarding the usage of `prima.py` could
be appropriate for the [Wien2k mailing list](http://www.wien2k.at/reg_user/mailing_list/).


CONTRIBUTING
------------

Bug reports and fixes, suggestions, and so on should be submitted via
GitHub as [Issues](https://github.com/eassmann/prima.py/issues) or
[Pull Requests](https://github.com/eassmann/prima.py/pulls).


ACKNOWLEDGEMENTS
----------------

Thanks to Maurits ⟨Maurits.Haverkort@cpfs.mpg.de⟩ for letting me build
upon his code.

Hex colors and color names are handled by the `webcolors` module,
written by James Bennett.

Thanks to Peter Blaha and Manuel Zingl for valuable comments on
previous versions.

`prima.py` was developed at [Vienna](http://www.ifp.tuwien.ac.at/cms)
and [Graz](https://itp.tugraz.at/) Universities of Technology.
