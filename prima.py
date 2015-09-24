#!/usr/bin/env python
# -*- coding: utf-8 -*-

## prima.py -- generate band-character plot using Maurits Haverkort's
##             SpaghettiPrimavera
##
## Copyright 2013 Elias Assmann <elias.assmann@gmail.com>

## prima.py version 0.1
##
## $Id: prima.py 76 2013-05-31 12:05:51Z assmann $

##
##	Vom Eise befreit sind Strom und Bäche
##	Durch des Frühlings holden, belebenden Blick,
##	Im Tale grünet Hoffnungsglück;
##	Der alte Winter, in seiner Schwäche,
##	Zog sich in rauhe Berge zurück.
##

from sppv import spaghettiprimavera, sppv_data as sppv

import os
import sys
import re

from getopt      import gnu_getopt
from numpy       import array,cumsum
from collections import deque

def svn_rev():
    rev = "$Rev: 76 $"

    try:
        return '-r' + re.search('(\d+)', rev).group()
    except:
        return ''

__version__ = "0.1" + svn_rev()

### The command line will be printed as a comment in the PS ###
sppv.cmdline[:] = ' '.join(sys.argv)

### Declare option related variables/functions ###
spin=''

def qtlname(x):
    global qtlfile
    qtlfile         = x
    sppv.qtlname[:] = x
def bandname(x):  sppv.bandname[:] = x
def efermi(x):    sppv.efermi      = x
def emin(x):      sppv.emin        = x
def emax(x):      sppv.emax        = x
def majortics(x): sppv.majortics   = x
def minortics(x): sppv.minortics   = x
def xsize(x):     sppv.xsize       = x
def ysize(x):     sppv.ysize       = x
def textsize(x):  sppv.textsize    = x
def fontname(x):  sppv.fontname[:] = x
def legend(x):    sppv.writelegend = True

thk0 = .1
def base_thk(x):
    global thk0
    thk0 = float(x)

def print_version(x):
    print 'prima.py version ' + __version__

    exit()

def print_help(x):
    print "USAGE: prima.py [OPTIONS] [IATM:IORB:R,G,B:THCK:LEG ...]"
    print "       (generate plot)"
    print "   OR: prima.py [-q QTL_BAND] [-s STRUCT]"
    print "       (print available atom and orbital names)"
    print
    print "   IATM: atom number(s) or name(s) N[+M...]"
    print "   IORB: orbital character number(s) or name(s) N[+M...]"
    print "         default: tot"
    print "  R,G,B: color"
    print "         default: 0,0,0 (black)"
    print "   THCK: line thickness increment"
    print "         default: 0 (no ``fat bands´´)"
    print "    LEG: legend entry"
    print "         default: IATM-IORB (if -L option given)"
    print
    print "OPTIONS:"

    for o in prima_options:
        print "\t-" + o[0] + ", --" + o[1] + "\t" + o[3]
    for o in sppv_options:
        print "\t-" + o[0] + ", --" + o[1] + "\t" + o[3]

    exit()

def do_up(x):
    global spin
    spin='up'

def do_dn(x):
    global spin
    spin='dn'

outfile=''
def outname(x):
    global outfile
    outfile = x

structfile = ""
def structname(x):
    global structfile
    structfile = x

sppv_options = [
    ('q:', 'qtl-file=',        qtlname,         'default: $(basename $PWD).qtl_band'),
    ('k:', 'klist-file=',      bandname,        'default: $(basename $PWD).klist_band'),
    ('s:', 'struct-file=',     structname,      'default: $(basename $PWD).struct     [optional]'),
    ('F:', 'fermi-energy=',    efermi,          'Fermi energy in Rydberg; default: from scf2'),
    ('e:', 'emin=',            emin,            '\tmin plotted energy in eV; default: -5'),
    ('E:', 'emax=',            emax,            '\tmax plotted energy in eV; default: +5'),
    ('t:', 'minor-tics=',      minortics,       '# minor tics per major interval; default: 4'),
    ('T:', 'major-tics=',      majortics,       'major tic distance in eV; default: 1'),
    ('X:', 'x-size=',          xsize,           '\tplot (not bbox) width  [pt]; default: 500'),
    ('Y:', 'y-size=',          ysize,           '\tplot (not bbox) height [pt]; default: 700'),
    ('S:', 'font-size=',       textsize,        'font size [pt?]; default: 12'),
    ('O:', 'font-name=',       fontname,        'font name; default: Times-Roman; Greek is in Symbol'),
    ('D:', 'base-thickness=',  base_thk,        'base line thickness [pt?]; default: 0.1'),
    ('L',  'write-legend',     legend,          'write legend'),
    ]

prima_options = [
    ('u',  'up',               do_up,            '\tsp: up'),
    ('d',  'dn',               do_dn,            '\tsp: down'),
    ('o:', 'out-file=',        outname,          'send output here instead of STDOUT'),
    ('h',  'help',             print_help,       '\tthis message'),
    ('v',  'version',          print_version,    'version info'),
    ]

### Parse options ###
opt_dict = {}
shopts = ""
longopts = []

for o in sppv_options + prima_options:
    l = o[1]
    if l[-1] == '=': l = l[0:-1]

    opt_dict['-' + o[0][0]] = o
    opt_dict['--' + l]      = o

    shopts += o[0]
    longopts += [o[1]]

(opts, args) = gnu_getopt(sys.argv[1:], shopts, longopts)

for (o, a) in opts:
    if opt_dict.has_key(o):
        opt_dict[o][2](a)

### Set defaults ###
case = os.getcwd().split(os.sep)[-1]

## set efermi from scf2 unless it has been set
if not [x for x in opts if x[0] in ['-F', '--fermi-energy']]:
    scf2 = case + '.scf2' + spin
    found = False

    if os.path.exists(scf2) and os.stat(scf2).st_size:
        for line in open(scf2, 'r'):
            if line.find(':FER') != -1:
                sppv.efermi = float(line.split()[-1])
                found = True
                break

    if not found:
        raise Exception("Fermi energy neither given nor found in " + scf2)

## set defaults for qtlfile, bandname, struct unless set
if not [x for x in opts if x[0] in ['-q', '--qtl-file']]:
    qtlfile          = case + '.qtl' + spin + '_band'
    sppv.qtlname[:]  = qtlfile
if not [x for x in opts if x[0] in ['-k', '--klist-file']]:
    sppv.bandname[:] = case + '.klist_band'
if not structfile:
    structfile = case + '.struct'


### Read QTL file for orbital character info ###
qtl = open(qtlfile, 'r')
jatomre = re.compile(' *JATOM *.* (tot[\d\w,]+)')
endre   = re.compile(' *BAND')

have_orbs = []

for line in qtl:
    m = jatomre.match(line)
    if m:
        o = {}

        for i,name in enumerate(m.group(1).split(',')):
            o[name] = i

        have_orbs.append(o)

    if endre.match(line): break

qtl.close()

have_orbs.append({'tot' : 0})   # interstitial


### Read struct file for atom info ###
# For each key (atom name, symbol, atomic number), have_atoms will
# contain a set of atom indices that correspond to that key
have_atoms = { 'all': set() }

try:
    struct = open(structfile, 'r')

    atomre = re.compile('(.*?) *NPT=')
    symre  = re.compile('[a-zA-Z]+')
    zre    = re.compile('Z: *([\d.]+)')

    iat = 1
    for line in struct:
        m = atomre.match(line)
        if m:
            # populate have_atoms: atom name
            aname = m.group(1)
            try:             have_atoms[aname].add(iat)
            except KeyError: have_atoms[aname] =  {iat}

            # atom symbol
            m = symre.match(aname)
            if m is None: continue

            asym = m.group()
            try:             have_atoms[asym].add(iat)
            except KeyError: have_atoms[asym] =  {iat}

            # atomic number
            m = zre.search(line)
            if m is None: continue

            az = "Z%g" % float(m.group(1))
            try:             have_atoms[az].add(iat)
            except KeyError: have_atoms[az] =  {iat}

            # 'all'
            have_atoms['all'].add(iat)

            iat += 1
    struct.close()
except IOError: pass            # this could be more robust

have_atoms['istl'] = {iat}


### Parse orbital character specs ###
atms = []
orbs = []
clrs = []
thks = []

legnames = []
legcolor = []

for a in args:
    tok = a.split(':')

    # legend entry provided, activate legend
    if len(tok) >= 5: legend(True)

    # defaults
    if len(tok) < 2: tok.append('tot')
    if len(tok) < 3: tok.append('0,0,0')
    if len(tok) < 4: tok.append(0)
    if len(tok) < 5: tok.append(tok[0]+'-'+tok[1])

    (aa, oo, c, t, l) = tok

    # transform color arg to array
    c = array(map(float, c.split(',')))

    # store legend entry -- there should be one entry for each cmd
    # line arg, not for each expanded atom/orb
    legnames.append(l)
    legcolor.append(c)

    # now expand atom / orbital args
    aa = deque(aa.split('+'))

    while aa:
        a = aa.popleft()

        try: iat = int(a)-1
        except ValueError:
            aa.extendleft(have_atoms[a])
            continue

        for o in oo.split('+'):
            try:               iorb = int(o)-1
            except ValueError: iorb = have_orbs[iat][o]
    
            atms.append(iat)
            orbs.append(iorb)
            clrs.append(c)
            thks.append(float(t))

if sppv.writelegend:
    # propagate legend info to Fortran
    sppv.colors = legcolor
    sppv.legend = ''.join(legnames)

    sppv.legentries = cumsum(map(len, legnames))

## no arguments: print available characters
if not atms:
    print 'prima.py using', qtlfile, 'and', structfile
    print
    for a,o in enumerate(have_orbs, start=1):
        print a, ':',
        for io in sorted(o.iteritems(), cmp=lambda a,b: a[1]-b[1]):
            print repr(io[0]) + ':' + repr(io[1]+1) + ',',
        else: print
    print
    for atn,iat in sorted(have_atoms.iteritems()):
        print atn, '=>', ','.join(map(repr, iat))

    sys.exit()

idx = range(len(atms))
istl = len(atms)

def chr2clr(orbchr):
    myclr = array([0., 0., 0.])

    for i in idx:
        myclr += clrs[i] * orbchr[atms[i], orbs[i]]#/(1-orbchr[istl, 0])

    return myclr


def chr2thk(orbchr):
    mythk = thk0;

    for i in idx:
        mythk += thks[i] * orbchr[atms[i],orbs[i]]#/(1-orbchr[istl, 0])

    return mythk


### Real Work ###

## redirect STDOUT if necessary -- need to use low-level I/O to
## propagate change to Fortran
if outfile:
    os.close(1)
    os.umask(~0644)
    os.open(outfile, os.O_WRONLY|os.O_CREAT|os.O_TRUNC) # should open on 1

spaghettiprimavera(chr2clr, chr2thk)
