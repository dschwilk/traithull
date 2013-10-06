traithull
=========

Provides an interface to the Qhull program that allows ecologists to easily calculate the convex hull volume (CHV) metric of functional diversity and do tests against null models. For a description of the method see [Cornwell et al 2006]

Author: Dylan Schwilk

Background
----------

Traithull provides an interface to the Qhull program of Barber and
Huhdanpaa (www.geom.umn.edu/software/qhull), which implements the
Quickhull algorithm for calculating convex hulls.  Traithull allows
easy input of species trait and plot data and provides randomization
output for testing null models of community assembly.

Traithull requires the qconvex executable from the qhull package. By
default, Traithull expects this program to be in the path. You may
also specify the location of qconvex on the traithull command line
with the -q option. The Microsoft Windows zip package of traithull
contains a copy of the qconvex.exe executable for windows

Installation
------------

### Linux


Install QHull. On Debian-based Linux distributions (Debian, Ubuntu,
etc) this is the package "qhull-bin". So a simple "apt-get install
qhull-bin" will work. Then unpack traithull and drop to a command line
in the traithull dierctory. If you want to install traithull so that
you can run the script from anywhere in your directory tree, copy
traithull.py to a directory in your path (~/bin for one-user
installation).

### Microsoft Windows

You need a working copy of QHull, specifically, make sure the
qconvex.exe is in your path or working directory. If you download the
windows zip file, it comes with qconvex.exe in the traithull
directory.

If you are not using the precompiled Microsoft Windows version, first
you must install qhull and add it to your path.  Then simply unpack
traithull to a desired location.

### Mac

You need a working python installation and a working qhull or at least
the qconvex executable. If someone wants to write up specific
instructions that would be great.


Command-line usage
------------------

To obtain a summary of command line usage, enter the following at the
shellcommand prompt:

   python traithull.py -h


To run the sample plot and species data, copy the examples to the main 
traithull directory and type the following at the command prompt:

   python traithull.py -pplots.txt species.txt
   
To produce random draws from the species pool, type

   python traithull.py species.txt


You can of course leave the example input files in the examples
directory and refer to them on the command line by path, eg
examples/plots.txt (or examples\plots.txt on Microsoft windows).


Command Summary
---------------

Output of "python traithull.py -h"

Usage: traithull.py [options] [trait_file]

'''
Options:
  --version             show program's version number and exit
  -h, --help            show this help message and exit
  -p PLOTFILE, --plotfile=PLOTFILE
                        File containing species occurance by plot
  -q QPATH, --qhull=QPATH
                        path to qconvex executable
  -r REPLICATES, --randsample=REPLICATES
                        Number of random samples per richness
  -d, --distance        Output mean and variance of nearest-neighbor distances
  -a, --Aussie          Output Aussie (Walker) fun. div. index
  -i, --individual      Do each treat in matrix individually (1-dimensional
                        version)
  -t, --total           Output total species pool results
  -v, --verbose         Verbose output
'''   
   
References
----------

[Cornwell et al 2006]: Cornwell, W.K, D.W. Schwilk and D.D. Ackerly. 2006. A trait-based test for habitat filtering: convex hull volume. _Ecology_ 87: 1465--1471.
