Version 0.3.0 (2015-11-xx):

 * new options to combine spins in one plot
   * `--mix-spins` forces the energies to be equal for both spins and
     mixes characters together to draw a single set of lines
   * `--join-spins` draws separate lines for each spin in one plot

 * new options `--out-suffix`, `--config-file`, `--verbose`

 * changed delimiter for combined atoms/orbitals to `++|&` to avoid
   clash with orbital names that contain `+`

 * make atoms and orbitals case insensitive

 * colors may be also given in HTML-style hex notation or color names;
   and may be prefixed by a factor, for scaling (`F × COLOR`)

 * read Fermi energy from `case.scf` rather than `case.scf2`; do not
   require it unless needed


Version 0.2.0 (2014-09-01):

 * added support for `x qtl`-generated files

 * fix some bugs related to parsing of band character names

 * added option `--case`

 * added line to postscript for explicit white background
   (current versions of `evince` were drawing black bg)

 * added optimization options in Makefile (can be important even here!)


Version 0.1 (2013-05-31): initial release
