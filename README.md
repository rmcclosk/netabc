# netabc

netabc is a computer program for fitting contact network models to transmission
trees by kernel approximate Bayesian computation.

## Prerequisites

You need the following shared libraries. 

* [GSL](http://www.gnu.org/software/gsl/)
* [Judy](http://judy.sourceforge.net/)
* [libyaml](http://pyyaml.org/wiki/LibYAML)
* [check](http://check.sourceforge.net/) (optional, for running tests)

On Ubuntu, these are all available in the repos.

   sudo apt-get install libgsl0ldbl libgsl0-dev 
   sudo apt-get install libjudydebian1 libjudy-dev 
   sudo apt-get install libyaml-0-2 libyaml-dev
   sudo apt-get install check

Optionally, to build the documentation, you need
[Doxygen](http://www.stack.nl/~dimitri/doxygen/index.html) and a LaTeX
distribution like TeX Live.

## Installation

If you just want to use the program, download one of the releases. The usual
procedure (./configure && make && make install) should work.

If you want to build straight from the git repo, for example if you are
interested in development, you need to run ./autogen.sh first to create
configure, then follow the same build procedure. This requires recent versions
of the autotools, flex, and bison.

## Use

### Available distributions

These are the distributions currently supported for specifying priors.
Currently, joint priors are not supported. In the case of _N_ and _I_ (the
total number of nodes, and the number of infected nodes), the priors you
specify will be truncated to the region _I_ ≤ _N_, and normalized so that each
marginal still integrates to 1. The "delta" distribution below is the Dirac
delta, equivalent to specifying an exact value.

Some distributions, such as the exponential distribution, are lower-bounded by
zero. The parameter called "shift" allows you to specify a different lower
bound, shifting the entire distribution. In other words, if the pdf of a
distribution is given by f(x), the shift parameter will modify the pdf to be
f(x) + shift. Similarly, the beta distribution is defined on the region [0, 1].
For this distribution, we introduce a parameter "scale" which scales the whole
distribution, so that it is defined on $[shift, shift + scale]$ instead. 

Since some distributions can be parameterized multiple ways, I tried to keep
the below table consistent with the Wikipedia page.

distribution | parameters
-------------|-----------
uniform | _a_ (lower limit), _b_ (upper limit)
gaussian | _μ_ (mean), _σ_ (variance)
delta | shift
exponential | shift, _λ_ (rate)
laplace | _μ_ (location), _b_ (scale)
exponential_power | _μ_ (location), α (scale), β (shape)
cauchy | _x<sub>0</sub>_ (location), _γ_ (scale)
rayleigh | shift, _σ_ (scale)
gamma | _k_ (shape), _θ_ (scale)
lognormal | _μ_ (location), _σ_ (scale) 
chi_squared | shift, _k_ (degrees of freedom)
f | shift, d<sub>1</sub> (first degrees of freedom), d<sub>2</sub> (second degrees of freedom)
student_t | shift, _ν_ (degrees of freedom)
beta | shift,  _α_ (first shape), _β_ (second shape), scale
logistic |  _μ_ (location), _s_ (scale)
pareto | _x<sub>m</sub>_ (scale), _α_ (shape)
weibull | shift, _λ_ (scale), _k_ (shape)
