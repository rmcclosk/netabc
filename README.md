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
interested in adding new network models, first clone the repository using
`git clone --recursive`. You must use `--recursive` because a forked version of
the igraph library has been packaged as a submodule. Then run `./autogen.sh` to
create the configure script, and build as usual. This requires recent versions
of the autotools, flex, and bison.

## Use

Netabc requires two inputs: an estimated transmission tree, in [Newick
format](https://en.wikipedia.org/wiki/Newick_format), and a
[YAML](http://yaml.org/) file specifying your priors on the network parameters.
The typical usage of netabc is as follows.
    
    netabc -d [output file] [other options] [tree] [yaml_file]

You can type `netabc -h` to get a summary of all the command line options. They
are described in detail in the next section. The structure of the YAML file
with priors is described in the following section.

Three additional binaries are included with `netabc`. The first is
`treekernel`, which computes the phylogenetic kernel of a pair of trees. The
second is `nettree`, which simulates a phylogeny over a transmission tree in
[GML](https://en.wikipedia.org/wiki/Graph_Modelling_Language) format. The third
is `treestat`, which computes several summary statistics on trees. Each of
these has a `--help` option which displays their usage.

### Command line options

**-t, --num-threads**

Number of threads to use. Ideally you would set this to the number of cores you
have available. The default is one thread, but we strongly encourage you to use
multiple threads if possible.

**-s, --seed**

An integer to use to set the state of the random number generator. Repeated
runs with the same tree and seed should produce the same results. The default
is to use the current time.

**-l, --decay-factor**

The decay factor parameter to the tree kernel (see [@poon2013mapping]). The
default is 0.3, and typical values are between 0.2 and 0.5. We encourage you to
perform a preliminary analysis to find optimal kernel meta-parameters for your
particular model, such as the cross-validation analysis performed in
[@poon2015phylodynamic].

**-g, --rbf-variance**

The variance for the radial basis function used in the tree kernel (see
[@poon2013mapping]). The default is 4, and typical values are between 0.5 and
8. As with the decay factor parameter, we encourage you to perform a separate
analysis to find an optimal value for your problem.

**-c, --nltt**

If this option is passed, multiply the tree kernel by the normalized
lineages-through-time (nLTT) statistic [@janzen2015approximate]. This is not
done by default.

**-n, --num-particles**

The number of particles to use for SMC (see [@del2012adaptive]). More particles
lead to a better SMC approximation, but at the expense of a linear increase in
computational cost. The default is 1000, which was used by [@del2012adaptive]
but leads to a fairly course approximation.

**-p, --num-samples**

The number of simulated trees to keep track of per particle. Again, more is
better in terms of accuracy, but the complexity is linearly related. The
default is 5; typical values are between 1 and 100.

**-q, --quality**

The coefficient of the equation which is solved to determine the next tolerance
value _ε_. It is called _α_ in [@del2012adaptive]. Essentially it represents a
tradeoff between speed and accuracy - lower values will finish with fewer
iterations but produce worse approximations. Typical values are between 0.9 and
0.99. The default is 0.95.

**-d, --trace**

Dump the values of each particle, their weights, and the distances of their
associated simulated trees to this file after each iteration. You *must*
specify this parameter if you want to record the estimated parameter values.

**-m, --net-type**

The name of the network model you want to fit. Currently available values are
listed in the next section. The default is "pa", which indicates
Barabasi-Albert preferential attachment.

**-e, --final-epsilon**

The first of two possible stopping conditions for SMC. The algorithm will be
stopped when the current distance tolerance falls to less than this value. The
default is 0.01. If you do not want to use this stopping criterion, pass `-e 0`
and specify a value for `-a`.

**-a, --final-accept-rate**

The second of two possible stopping conditions for SMC. The algorithm will be
stopped when the acceptance proportion of MCMC moves of the particles falls below
this value. The default is 0.015, or 1.5%. If you do not wish to use this
stopping criterion, pass `-a 0` and specify a value for `-e`.

### Specifying the model and priors

Currently, three models are supported.

model | description
------|------------
pa    | Barabasi-Albert preferential attachment (uses `igraph_barabasi_game`)
gnp   | Erdos-Renyi random network (uses `igraph_erdos_renyi_game`)
smallworld   | Watts-Strogatz small world network (uses `igraph_watts_strogatz_game`)

Priors are specified in a YAML file. It's also possible to specify exact values
on the parameters. Consider the following example YAML file.

    N: ["uniform", 287, 10000]
    I: ["uniform", 287, 10000]
    time: 0
    transmit_rate: 1
    remove_rate: 0
    m: ["uniform", 2, 6]
    alpha: ["uniform", 0, 2]

Here, we are specifying uniform priors with various endpoints for the
parameters "N", "I", "m", and "alpha". The parameters "time", "transmit_rate",
and "remove_rate" have been specified exactly at the values 0, 1, and 0,
respectively.

You must specify either priors for all parameters of the model you are fitting.

model | parameter | meaning
------|-----------|--------
all   | N         | number of nodes in the network
all   | I         | number of infected nodes at time of sampling
all   | time      | amount of simulated time passed at time of sampling
all   | transmit_rate | transmission rate per discordant edge
all   | remove_rate | removal rate per infected node
pa    | m         | number of edges added per vertex
pa    | alpha     | preferential attachment power
gnp   | p         | probability per edge
smallworld | nei  | number of neighbours for each vertex
smallworld | p    | rewiring probability

The required parameters for all models have some dependency on each other. In
the case of _N_ and _I_ (the total number of nodes, and the number of infected
nodes), the priors you specify will be truncated to the region _I_ ≤
_N_ by rejection sampling. If there is not enough prior density in the region 
_I_ ≤ _N_, the initialization of the particles will fail.

The time parameter indicates the simulation time at which sampling occurs. If
set to 0, it is ignored, and the simulation proceeds as long as necessary _I_
nodes to become infected. Likewise, if _I_ is set to 0, the simulation proceeds
until the specified time has been reached regardless of how large the epidemic
gets. If _I_ and time are both 0, then the simulation will always continue
until there are no more discordant edges in the network. You cannot specify
exact values for time or _I_ unless the other is zero.

These are the distributions currently supported for specifying priors.
Currently, joint priors are not supported, so you must list priors for each
parameter independently. The "delta" distribution is the Dirac delta,
equivalent to specifying an exact value.

Some distributions, such as the exponential distribution, are lower-bounded by
zero. The parameter called "shift" allows you to specify a different lower
bound, shifting the entire distribution. In other words, if the pdf of a
distribution is given by f(x), the shift parameter will modify the pdf to be
f(x) + shift. Similarly, the beta distribution is defined on the region [0, 1].
We introduce a parameter "scale" which scales the whole distribution, so that
it is defined on $[shift, shift + scale]$ instead. 

Some distributions can be parameterized multiple ways. I tried to keep the
below table consistent Wikipedia.

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

## Adding new models

To maintain the speed of the program, it is necessary to directly modify the C
code if you want to add a new network model. All the model-specific code is in
`src/netabc.c`. At this stage, the procedure to add a new model is somewhat
involved.

1. Add a new element to the `net_type` enum for your new model.
2. Append the number of model parameters it has to the global `NUM_PARAMS`
   array.
3. Add elements to the `net_parameter` enum for all the model's parameters.
4. Add the parameters' names to the `PARAM_NAMES` array.
5. Modify the `get_options()` function to accept an abbreviation for your
   network type on the command line. Look for the line `opts.net = NET_TYPE_PA`
   and add something similar. Add the abbreviation to the `usage()` function.
6. Add code to simulate a network under your model in the `sample_dataset()`
   function. Look for the line `case NET_TYPE_PA:` and add something similar.
