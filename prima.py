#!/usr/bin/env python
# -*- coding: utf-8 -*-

## prima.py -- generate band-character plot using Maurits Haverkort's
##             SpaghettiPrimavera
##
## Copyright 2013-2015 Elias Assmann <elias.assmann@gmail.com>

##
##	Vom Eise befreit sind Strom und Bäche
##	Durch des Frühlings holden, belebenden Blick,
##	Im Tale grünet Hoffnungsglück;
##	Der alte Winter, in seiner Schwäche,
##	Zog sich in rauhe Berge zurück.
##

from sppv import spaghettiprimavera, sppv_data as sppv

import os, sys, re, getopt, numpy, collections, webcolors, contextlib, \
  traceback, subprocess, warnings

def git_version():
    rev = "$version:v0.2.0-29-g582dece$"

    try:
        return re.search(":\s*(.*?)\s*\$", rev).group(1)
    except:
        return 'unknown'


__version__ = git_version()
sppv.version[:] = __version__

def carp(thing):
    if verbose and isinstance(thing, Exception):
        traceback.print_exc()
    else:
        if isinstance(thing, tuple):
            thing = ' '.join(str(x) for x in thing)
        print >>sys.stderr, "prima.py:", thing

def croak(thing, status=1):
    carp(thing)
    sys.exit(status)

@contextlib.contextmanager
def open_w_error(filename, mode="r"):
    try:
        f = open(filename, mode)
    except IOError, err:
        yield None, err
    else:
        try:
            yield f, None
        finally:
            f.close()

def find_reverse(lefi, string):
    """Find first line containing ‘string’ in ‘lefi’, searching backwards.
    ‘lefi’ must be seekable."""

    pos = lefi.tell()
    line  = ''
    while pos > 0:
        pos -= 1
        lefi.seek(pos)

        c = lefi.read(1)
        if c=='\n':
            if string in line:
                return line

            line = ''
        else:
            line = c + line

    if string in line: return line
    else:              return None

### The command line will be printed as a comment in the PS ###
sppv.cmdline[:] = ' '.join(sys.argv)

### Declare option related variables/functions ###
bandfile = None
def bandname(x):
    global bandfile
    print>>log, 'set klist file = ‘'       +x+'’'
    bandfile            = x
    sppv.bandname[:]    = x
def efermi(x):
    print>>log, 'set Fermi energy = ‘'     +x+'’ Ryd'
    sppv.efermi         = x
def emin(x):
    print>>log, 'set bottom energy = ‘'    +x+'’ eV'
    sppv.emin           = x
def emax(x):
    print>>log, 'set top energy = ‘'       +x+'’ eV'
    sppv.emax           = x
def majortics(x):
    print>>log, 'set major tics every = ‘' +x+'’ eV'
    sppv.majorticks     = x
def minortics(x):
    print>>log, 'set minor tics = ‘'       +x+'’ per major tic interval'
    sppv.minorticks     = x
def xsize(x):
    print>>log, 'set x-size = ‘'           +x+'’ pt'
    sppv.xsize          = x
def ysize(x):
    print>>log, 'set y-size = ‘'           +x+'’ pt'
    sppv.ysize          = x
def textsize(x):
    print>>log, 'set font size = ‘'        +x+'’ pt'
    sppv.textsize       = x
def fontname(x):
    print>>log, 'set font = ‘'             +x+'’'
    sppv.fontname[:]    = x
def legend(x):
    print>>log, 'set write legend'
    sppv.writelegend    = True
def nolegend(x):
    print>>log, 'unset legend'
    sppv.writelegend    = False
def noklines(x):
    print>>log, 'unset k-lines'
    sppv.writeklines    = False
def klines(x):
    print>>log, 'set draw k-lines'
    sppv.writeklines    = True
def noklabels(x):
    print>>log, 'unset k-labels'
    sppv.writeklabels   = False
def klabels(x):
    print>>log, 'set write k-labels'
    sppv.writeklabels   = True
def nofermi(x):
    print>>log, 'unset Fermi energy line'
    sppv.writefermiline = False
def fermi(x):
    print>>log, 'set draw Fermi energy line'
    sppv.writefermiline = True
def axes_thk(x):
    print>>log, 'set axis thickness = ‘'+x+'’ pt'
    sppv.axesthickness  = x
def base_thk(x):
    print>>log, 'set line base thickness = ‘'+x+'’ pt'
    sppv.basethickness  = x

qtlfile  = None
qtlfile1 = None
def qtlname(x):
    global qtlfile, qtlfile1
    qtlfile1        = qtlfile
    qtlfile         = x
    sppv.qtlname[:] = x
    print>>log, 'set qtl file = ‘'+x+'’'

def print_version(x):
    print 'prima.py version ' + __version__

    sys.exit()

def print_help(x):
    print "USAGE: prima.py [OPTIONS] [IATM:IORB:R,G,B:THCK:LEG ...]"
    print "       (generate plot)"
    print "   OR: prima.py [-q QTL_BAND] [-s STRUCT]"
    print "       (print available atom and orbital names)"
    print
    print "   IATM: atom number(s) or name(s) N[&M...]"
    print "   IORB: orbital character number(s) or name(s) N[&M...]"
    print "         default: tot"
    print "  R,G,B: color"
    print "         default: 0,0,0 (black)"
    print "   THCK: line thickness increment"
    print "         default: 0 (no “fat bands”)"
    print "    LEG: legend entry"
    print "         default: IATM-IORB (if -L option given)"
    print
    print "OPTIONS:"

    for o in [o for o in prima_options if o[3] is not None]:
        print "\t-" + o[0] + ", --" + o[1] + "\t" + o[3]
    for o in [o for o in sppv_options  if o[3] is not None]:
        print "\t-" + o[0] + ", --" + o[1] + "\t" + o[3]
    print
    print "Boolean options as listed negate the defaults; long options to reinstate"
    print "the defaults (e.g. ‘--no-legend’, ‘--k-labels’) are also valid."

    sys.exit()

spin = ''
mix  = False
def do_up(x):
    global spin
    spin='up'
    print>>log, 'set spin = ‘'+spin+'’'
def do_dn(x):
    global spin
    spin='dn'
    print>>log, 'set spin = ‘'+spin+'’'
def no_spin(x):
    global spin
    spin=''
def do_mix(x):
    global spin, mix
    spin = 'updn'
    mix  = True
    print>>log, 'mix spins'
    print>>log, 'set spin = ‘'+spin+'’'
def do_join(x):
    global spin, mix
    spin = 'updn'
    mix  = False
    print>>log, 'join spins'
    print>>log, 'set spin = ‘'+spin+'’'

class nofile(object):
    def write(self, string):
        pass

verbose = False
log = nofile()
def do_verbose(x):
    global verbose, log
    verbose = True
    log = sys.stderr
    print>>log, 'set verbose'
def no_verbose(x):
    global verbose, log
    verbose = False
    print>>log, 'unset verbose'
    log = nofile()

outsuffix = None
def suffixname(x):
    global outsuffix
    outsuffix = x
    print>>log, 'set output suffix = ‘'+x+'’'
    print>>log

outfile = None
def outname(x):
    global outfile, outsuffix
    outfile = x
    outsuffix = None           # delete outsuffix for explicit outfile
    print>>log, 'set output file = ‘'+x+'’'

structfile = None
def structname(x):
    global structfile
    structfile = x
    print>>log, 'set struct = ‘'+x+'’'

scffile = None
def scfname(x):
    global scffile
    scffile = x
    print>>log, 'set scf file = ‘'+x+'’'

case = os.getcwd().split(os.sep)[-1]
print>>log, case
def casename(x):
    global case
    case = x
    print>>log, 'set case = ‘'+x+'’'

def confname(x): # this is a no-op because config files are special
    # but it is a convenient place to verbosely echo config files
    print>>log, 'config file ‘'+x+'’read'
    

sppv_options = [
    ('f:', 'case=',           casename,     '\tdefault: $(basename $PWD)'),
    ('q:', 'qtl-file=',       qtlname,      'default: ‘case’.qtl[sp][_band]'),
    ('k:', 'klist-file=',     bandname,     'default: ‘case’.klist_band'),
    ('s:', 'struct-file=',    structname,   'default: ‘case’.struct\t\t[optional]'),
    ('c:', 'scf-file=',       scfname,      'default: ‘case’.scf\t\t[optional, for Fermi enerfy]'),
    ('F:', 'fermi-energy=',   efermi,       'Fermi energy in Rydberg; default: from scf'),
    ('e:', 'emin=',           emin,         '\tmin plotted energy in eV; default: -5'),
    ('E:', 'emax=',           emax,         '\tmax plotted energy in eV; default: +5'),
    ('t:', 'minor-tics=',     minortics,    '# minor tics per major interval; default: 4'),
    ('T:', 'major-tics=',     majortics,    'major tic distance in eV; default: 1'),
    ('X:', 'x-size=',         xsize,        '\tplot (not bbox) width  [pt]; default: 500'),
    ('Y:', 'y-size=',         ysize,        '\tplot (not bbox) height [pt]; default: 700'),
    ('S:', 'font-size=',      textsize,     'font size [pt?]; default: 12'),
    ('O:', 'font-name=',      fontname,     'font name; default: Times-Roman; Greek is in Symbol'),
    ('D:', 'base-thickness=', base_thk,     'base line thickness [pt?]; default: 0.1'),
    ('A:', 'axes-thickness=', axes_thk,     'line thickness for axes [pt?]; default: 0.5'),
    ('L',  'legend',          legend,       '\twrite legend'),
    (None, 'no-legend',       nolegend,     None),
    ('K',  'no-k-labels',     noklabels,    'omit labels for k-points'),
    (None, 'k-labels',        klabels,      None),
    ('G',  'no-k-lines',      noklines,     'omit vertical lines for labeled k-points'),
    (None, 'k-lines',         klines,       None),
    ('0',  'no-fermi-line',   nofermi,      'omit horizontal line for Fermi energy'),
    (None, 'fermi-line',      fermi,        None),
]

prima_options = [
    ('u',  'up',              do_up,         '\tsp: up'),
    ('d',  'dn',              do_dn,         '\tsp: down'),
    ('m',  'mix-spins',       do_mix,        '\tmix spins  [equal energies forced]'),
    ('j',  'join-spins',      do_join,       'join spins [energies may be different]'),
    ('N',  'no-spin',         no_spin,       '\tnon-spin-polarized'),
    ('o:', 'out-file=',       outname,       'send output here instead of STDOUT'),
    ('a:', 'out-suffix=',     suffixname,    'send output to $case$suffix instead of STDOUT'),
    ('C:', 'config-file=',    confname,      'read options from here'),
    ('h',  'help',            print_help,    '\tthis message'),
    ('v',  'verbose',         do_verbose,    '\tprint additional info to STDERR'),
    (None, 'no-verbose',      do_verbose,    None),
    ('V',  'version',         print_version, '\tversion info'),
]

### Parse options ###
opt_dict = {}
shopts = ""
longopts = []

for o in sppv_options + prima_options:
    l = o[1]
    if l[-1] == '=': l = l[0:-1]

    opt_dict[l] = o
    longopts += [o[1]]

    if o[0] is not None:
        opt_dict[o[0][0]] = o
        shopts   +=  o[0]

print>>log, 'Parsing option arguments …'

try:
    (opts, args) = getopt.gnu_getopt(sys.argv[1:], shopts, longopts)
except getopt.GetoptError as e:
    croak(e)

## First, look for and read config files
xopts = []
skipre = re.compile('\s*(#.*)?$')
for (o, a) in opts:
    if o not in ['-C', '--config-file']: continue

    try:
        c = open(a, 'U')

        for line in c:
            if skipre.match(line): continue

            s = line.rstrip('\n').split(None, 1)

            lo = s[0]

            if len(s) > 1:
                la = s[1]
            else:
                la = ''

            xopts += [(lo, la)]

    except IOError as e:
        croak(('error reading config file:', e))

opts = xopts + opts
    
## Now for the other options
stripdash = re.compile('^--?')

for (o, a) in opts:
    o = stripdash.sub('', o)
    opt_dict[o][2](a)

### Set defaults ###

## This is done here because it relies on ‘case’ and some other
## options being set

if args:
    ## set efermi from scf unless it has been set
    if not [x for x in opts if x[0] in ['-F', '--fermi-energy']]:
        try:
            if scffile is None: scffile = case + '.scf'

            print>>log, 'no -F option, will try to find Fermi energy in ‘'\
                         +scffile+'’'

            with open(scffile, 'U') as scf:
                scf.seek(0, os.SEEK_END)

                line = find_reverse(scf, ':FER')
                if spin=='up':
                    line = find_reverse(scf, ':FER')

                if line is None:
                    raise Exception('Fermi energy not found in ‘' +
                                    scffile + '’')

                sppv.efermi = float(line.split()[-1])

        except Exception as e:
            croak(('error trying to fetch Fermi energy:', e))

    print>>log, 'Fermi energy is EF =', sppv.efermi, 'Ryd'

## set defaults for qtlfile, bandname, struct unless set
print>>log
print>>log, 'Files used:'
if qtlfile    is None:
    s = spin
    t = 'up'
    if spin == 'updn': s = 'dn'

    qtlfile      = case + '.qtl' + s + '_band'
    qtlfile1     = case + '.qtl' + t + '_band'
    if not os.path.exists(qtlfile):
        qtlfile  = case + '.qtl' + s
        qtlfile1 = case + '.qtl' + t

    sppv.qtlname[:] = qtlfile

    if spin == 'updn':
        print>>log, '   qtl1= ‘'+qtlfile1+'’'
    print>>log,     '   qtl = ‘'+qtlfile +'’'
elif spin=='updn' and  qtlfile1 is None:
    croak('need 0 or 2 qtl files in updn mode, not 1')

if bandfile   is None:
    bandfile         = case + '.klist_band'
    sppv.bandname[:] = bandfile
    print>>log, '   klist_band = ‘'+bandfile+'’'
if structfile is None:
    structfile = case + '.struct'
    print>>log, '   struct = ‘'+structfile+'’'

## set outfile from outsuffix if set
if outsuffix is not None:
    outfile = case + outsuffix

if outfile is not None:
    print>>log, '   output to ‘'+outfile+'’'

### Read QTL file for orbital character info ###
print>>log
print>>log, 'Reading qtl file …'
try:
    qtl     = open(qtlfile, 'U')
except IOError as e:
    croak(('error opening qtl file:', e))

jatomre = re.compile(' *JATOM *.* (tot.*)')
endre   = re.compile(' *BAND')

have_orbs = []
for line in qtl:
    m = jatomre.match(line)
    if m:
        print>>log
        print>>log, 'Next atom'
        o = {}

        for i,name in enumerate(m.group(1).split(',')):
            name = name.strip()

            if name == '': break

            name = {
                '0': 'S',
                '1': 'P',
                '2': 'D',
                '3': 'F'
                }.get(name, name).lower()
            o[name] = i
            print>>log, '   found orbital', name, '=>', i

        have_orbs.append(o)

    if endre.match(line): break

qtl.close()

num_atoms = len(have_orbs)

have_orbs.append({'tot' : 0})   # interstitial


### Read struct file for atom info ###
# For each key (atom name, symbol, atomic number), have_atoms will
# contain a set of atom indices that correspond to that key
have_atoms = { 'all': set() }

try:
    print>>log
    print>>log, 'Reading struct file …'
    struct = open(structfile, 'U')

    atomre = re.compile('(.*?) *NPT=')
    symre  = re.compile('[a-zA-Z]+')
    zre    = re.compile('Z: *([\d.]+)')

    iat = 1
    for line in struct:
        m = atomre.match(line)
        if m:
            # populate have_atoms: atom name
            aname = m.group(1)
            try:             have_atoms[aname.lower()].add(iat)
            except KeyError: have_atoms[aname.lower()] =  {iat}

            # atom symbol
            m = symre.match(aname)
            if m is None: continue

            asym = m.group()
            try:             have_atoms[asym.lower()].add(iat)
            except KeyError: have_atoms[asym.lower()] =  {iat}

            # atomic number
            m = zre.search(line)
            if m is None: continue

            az = "z%g" % float(m.group(1))
            try:             have_atoms[az].add(iat)
            except KeyError: have_atoms[az] =  {iat}

            # 'all'
            have_atoms['all'].add(iat)

            print>>log, '   found atom', iat, '<=', aname, asym, az

            iat += 1

    have_atoms['istl'] = {iat}

    if iat-1 != num_atoms:
        croak('struct (%i atoms) and qtl (%i atoms) files inconsistent'
              %(iat-1, num_atoms))

    struct.close()
except IOError: pass                # FIXME: this could be more robust

### Parse orbital character specs ###
atms = []
orbs = []
clrs = []
thks = []

legnames  = []
legcolors = []

combinator = re.compile('\+\+|&')
multiplier = re.compile('x|×')

def parsecolor(cc):
    if multiplier.search(cc):
        (fact, cc) = multiplier.split(cc, 1)
        fact = float(fact)
    else:
        fact = 1.

    cc = cc.strip()

    try:
        c = numpy.array(map(float, cc.split(',')))

        if c.size == 1: c = c.repeat(3)
        if c.size != 3: raise StandardError()

        return c * fact
    except ValueError: pass

    try:
        c = webcolors.hex_to_rgb(cc)
        return numpy.array(c)/255. * fact
    except ValueError: pass

    c = webcolors.name_to_rgb(cc)
    return numpy.array(c)/255. * fact

print>>log
print>>log, 'Parsing non-option arguments …'

for a in args:
    tok = a.split(':', 4)

    # legend entry provided, activate legend
    if len(tok) > 4: legend(True)

    tok += [''] * (5 - len(tok))

    print>>log, '   found tokens', tok,

    (aa, oo, cc, t, l) = tok[:5]

    # defaults
    if oo == '': oo = 'tot'
    if cc == '': cc = '0,0,0'
    if t  == '': t  = '0'
    if l  == '': l  = tok[0]+'-'+oo
    
    aa = aa.lower()
    oo = oo.lower()

    # transform color arg to array
    try:
        c = parsecolor(cc)
    except StandardError:
        croak(('invalid color spec:', cc))

    # store legend entry -- there should be one entry for each command
    # line arg, not for each expanded atom/orbital
    if l:
        legnames .append(l)
        legcolors.append(c)

    print>>log, '=> [', aa, oo, c, t, l, ']'

    # now expand atom / orbital args
    aa = collections.deque(combinator.split(aa))

    upsym = ['-up', '-↑', '↑']
    dnsym = ['-dn', '-↓', '↓']
    udsym = ['-ud', '-↓↑', '↑↓']
    spinre = re.compile('(.*)' + '(' + '|'.join(upsym+dnsym+udsym) + ')')
    mult = 2 if spin=='updn' else 1
    while aa:
        a = aa.popleft()

        add = 0
        sp  = ''
        if spin=='updn' and isinstance(a, str):
            m = spinre.match(a)
            if m:
                (a, sp) = m.groups()
            else:                        # empty spin part means up+dn
                sp = udsym[0]

        if sp in dnsym: add = 1

        if sp in udsym:
            aa += [a+upsym[0], a+dnsym[0]]
            continue

        try:
            iat = int(a)-1
        except ValueError:
            a = a.strip()
            try:
                xa = [str(i)+sp for i in have_atoms[a]]
                aa.extendleft(xa)
                continue
            except KeyError:
                croak('atom ‘' + a + '’ not found')

        for o in combinator.split(oo):
            o = o.strip()

            try:
                iorb = int(o)-1

                if iorb < 0 or iorb >= len(have_orbs[iat]):
                    raise KeyError(iorb)
            except ValueError:
                try:
                    iorb = have_orbs[iat][o]
                except KeyError:
                    croak('orbital ‘' + o + '’ not found')
            except KeyError:
                croak('orbital ‘' + o + '’ not found')

            atms.append(1 + mult*iat + add)      # +1 for Fortran indices
            orbs.append(1 + iorb)
            clrs.append(c)
            thks.append(float(t))

            print>>log, \
                '      adding atom %i, orb %i, color %s, thickness %g' \
                %(1+2*iat+add, 1 + iorb, c, float(t))

print>>log
print>>log, "Atoms   :", atms
print>>log, "Orbitals:", orbs
print>>log, "Colors  :", clrs
print>>log, "Thickn's:", thks

sppv.orbcolors = numpy.array(clrs).T
sppv.atomsorbs = numpy.array((atms, orbs))
sppv.orbthicks = thks

if sppv.writelegend:
    # propagate legend info to Fortran
    sppv.legcolors = numpy.array(legcolors).T

    # FIXME: the ‘+ ' '’ seems to prevent a segfault if a single
    # legend entry of one character is given
    sppv.legend = ''.join(legnames) + ' '

    sppv.legentries = numpy.cumsum(map(len, legnames))

## no arguments: print available characters
if not atms:
    print>>log
    print>>log, 'No non-option arguments: print available characters'

    print 'prima.py', __version__, 'using', qtlfile, 'and', structfile
    print
    for a,o in enumerate(have_orbs, start=1):
        print a, ':',
        for io in sorted(o.iteritems(), cmp=lambda a,b: a[1]-b[1]):
            print repr(io[0]) + ':' + repr(io[1]+1) + ',',
        else: print
    print
    for atn,iat in sorted(have_atoms.iteritems()):
        atn = atn[0].upper() + atn[1:]
        print atn, '=>', ','.join(map(repr, iat))

    print>>log
    print>>log, 'All done.'
    sys.exit()

### Real Work ###

## In ‘updn’ mode, open a pipe to paste the two qtl files together.  
if spin=='updn':
    print>>log
    print>>log, 'Opening pipe to paste together ‘up’ and ‘dn’ qtl files …'

    warnings.simplefilter('ignore')
    qtlpipe = os.tempnam(None, 'udqtl')
    warnings.resetwarnings()
    try:
        os.mkfifo(qtlpipe)
        sppv.qtlname[:] = qtlpipe

        if os.fork()==0: # child
            # open() will block until fifo opened for reading
            with open(qtlpipe, 'w') as qtl, \
                 open(qtlfile1) as q1, open(qtlfile) as q2:
                # header
                natre = re.compile('NAT\s*=\s*(\d+)')
                l1=''
                while 'BAND' not in l1:
                    print>>qtl, l1,

                    if 'JATOM' in l1:
                        print>>qtl, l2,

                    m = natre.match(l1)
                    if m:
                        if int(m.group(1)) != num_atoms:
                            croak(('inconsistent number of atoms in ‘%s’ ' +
                                   '(%i atoms instead of %i)'
                                   %(qtlfile1, int(m.group(1)), num_atoms)))

                        m = natre.match(l2)
                        if int(m.group(1)) != num_atoms:
                            croak(('inconsistent number of atoms in ‘%s’ ' +
                                   '(%i atoms instead of %i)'
                                   %(qtlfile2, int(m.group(1)), num_atoms)))
                    # get new lines
                    l1 = q1.readline()
                    l2 = q2.readline()

                # header done, l1 and l2 contain “BAND 1”; print extra
                # line for the extra interstitial
                print>>qtl, 'JATOM X interstitial'

                nk=0
                while l1 and l2:
                    nk += 1
                    print>>qtl, l1,

                    l1 = q1.readline()

                    # handling of the q2 file depends on mix/join
                    if mix:
                        l2 = q2.readline()
                        if not l2.startswith(' BAND'):
                            print>>qtl, l2,
                    else:
                        if not l1 or l1.startswith(' BAND'):
                            print>>qtl, l2,
                            for k in range(nk-1):
                                l2 = q2.readline()
                                print>>qtl, l2,
                                print>>qtl, l2[0:14], ' 0'*20
                            nk=0
                            l2 = q2.readline()
                        else:
                            print>>qtl, l1[0:14], ' 0'*20


            os.unlink(qtlpipe)
            sys.exit()
    except EnvironmentError as e:
        croak(('error opening pipe for --updn qtl:', e))

    print>>log, 'Using named pipe ‘'+qtlpipe+'’as qtl file'

## Redirect STDOUT if necessary -- need to use low-level I/O to
## propagate change to Fortran
if outfile not in [None, '-']:
    print>>log
    print>>log, 'Redirecting STDOUT for Fortran'
    os.close(1)
    os.umask(~0644)
    # should open on 1 (fingers crossed)
    fid = os.open(outfile, os.O_WRONLY|os.O_CREAT|os.O_TRUNC)
    if fid != 1:
        croak(('Redirecting STDOUT did not work -- opened on', fid,
               'instead of 1'))

print>>log
print>>log, 'Handing over to Fortran …'

spaghettiprimavera()

print>>log
print>>log, 'All done.'


## Time-stamp: <2015-09-22 16:46:36 assman@faepop36.tu-graz.ac.at>
