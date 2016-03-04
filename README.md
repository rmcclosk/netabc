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

## Use

### Available distributions

These are the distributions currently supported for specifying priors.
Currently, joint priors are not supported. In the case of $N$ and $I$ (the
total number of nodes, and the number of infected nodes), the priors you
specify will be truncated to the region $I \leq N$, and normalized so that each
marginal still integrates to 1. The "delta" distribution below is the Dirac
delta, equivalent to specifying an exact value.

Some distributions, such as the exponential distribution, are canonically
lower-bounded by zero. The parameter called "shift" allows you to specify a
different lower bound, shifting the entire distribution. In other words, if the
pdf of a distribution is given by f(x), the shift parameter will modify the pdf
to be f(x) + shift.

Similarly, the beta distribution is defined on the region $[0, 1]$. For this
distribution, we introduce a parameter "scale" which scales the whole
distribution, so that it is defined on $[shift, shift + scale]$ instead. 

Since some distributions can be parameterized multiple ways, we tried to keep
the below table consistent with the Wikipedia page.

distribution | parameters
-------------|-----------
uniform | $a$ (lower limit), $b$ (upper limit)
gaussian | $\mu$ (mean), $\sigma$ variance
delta | shift
exponential | shift, $\lambda$ (rate)
laplace | $\mu$ (location), $b$ (scale)
exponential_power | $\mu$ (location), $\alpha$ (scale), $\beta$ (shape)
cauchy | $x_0$ (location), $\gamma$ (scale)
