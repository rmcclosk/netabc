# netabc

netabc is a computer program for fitting contact network models to transmission
trees by kernel approximate Bayesian computation.

## Prerequisites

You need the following shared libraries. 

* [GSL](http://www.gnu.org/software/gsl/)
* [Judy](http://judy.sourceforge.net/)
* [libyaml](http://pyyaml.org/wiki/LibYAML)
* [check](http://check.sourceforge.net/)
* [igraph (my fork)](https://github.com/rmcclosk/igraph)

You also need the pthreads API, and a fairly recent version of Flex and Bison.
If you're on unix these are probably already installed. If you're not on unix,
this project probably won't work for you :)

Optionally, to build all the documentation, you need
[Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) and a LaTeX
distribution like TeX Live.

## Installation

The usual procedure (./configure && make && make install) should work, although
I'd recommend making a separate build directory and installing from there. If
the repo didn't come with a configure script (I'll be adding it later when
things aren't in as much flux), do ./autogen.sh first to create it.
