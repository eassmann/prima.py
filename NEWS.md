Version 0.3 (2015-11-xx):

 * add options to combine spins in one plot
   * `--mix-spins` forces the energies to be equal for both spins and
     mixes characters together to draw a single set of lines
     (may be better for antiferromagnetism)
   * `--join-spins` draws separate lines for each spin in one plot
     (more appropriate for ferromagnetism)

 * read Fermi energy from `case.scf` rather than `case.scf2`; do not
   require it unless needed

 * change delimiter for combined atoms/orbitals to `++|&` to avoid
   clash with orbital names that contain `+`

 * new options `--out-suffix`, `--config-file`, `--verbose`

 * make atoms and orbitals case insensitive

 * colors may be also given in HTML-style hex notation or color names;
   and may be prefixed by a factor, for scaling (`F × COLOR`)


Version 0.2 (2014-09-01):

 * add support for `x qtl`-generated files

 * fix some bugs related to parsing of band character names

 * added option `--case`

 * add line to postscript for explicit white background
   (current versions of `evince` were drawing black bg)

 * add optimization options in Makefile (can be important even here!)


Version 0.1 (2013-05-31): initial release
