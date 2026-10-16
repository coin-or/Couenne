# Couenne

[![A COIN-OR Project](https://coin-or.github.io/coin-or-badge.png)](https://www.coin-or.org)

Projects such as this one are maintained by a small group of volunteers under
the auspices of the non-profit [COIN-OR Foundation](https://www.coin-or.org)
and we need your help! Please consider [sponsoring our
activities](https://github.com/sponsors/coin-or) or [volunteering](mailto:volunteer@coin-or.org) to help!

[![Latest Release](https://img.shields.io/github/v/release/coin-or/Couenne?sort=semver)](https://github.com/coin-or/Couenne/releases)

_This file is auto-generated from [config.yml](.coin-or/config.yml) using the 
[generate_readme](.coin-or/generate_readme) script.
To make changes, please edit [config.yml](.coin-or/config.yml) or the generation scripts
[here](.coin-or/generate_readme) and [here](https://github.com/coin-or/coinbrew/blob/master/scripts/generate_readme)._

Couenne (*C*onvex *O*ver and *U*nder *EN*velopes for *N*onlinear *E*stimation)
is a branch-and-bound algorithm to solve Mixed-Integer Nonlinear Programming
(MINLP) problems of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; min <i>f</i><sub>0</sub>(<i>x</i>,<i>y</i>)<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     <i>f</i><sub>i</sub>(<i>x</i>,<i>y</i>) &le; 0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>i</i>=1,2..., <i>m</i><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     <i>x</i> in R<sup>n</sup>,
                                   <i>y</i> in Z<sup>p</sup>

where all _f_<sub>i</sub>(_x_,_y_) are, in general, nonlinear functions.

Couenne aims at finding global optima of nonconvex MINLPs. It implements
linearization, bound reduction, and branching methods within a
branch-and-bound framework. Its main components are:
 * an expression library;
 * separation of linearization cuts;
 * branching rules;
 * bound tightening methods.

Couenne builds on top of [Bonmin](https://github.com/coin-or/Bonmin),
[Cbc](https://github.com/coin-or/Cbc), and [Ipopt](https://github.com/coin-or/Ipopt).
It can be used as a callable library or as a stand-alone executable that reads
models written in [AMPL](https://ampl.com).

Couenne can also be used through
 * AMPL (natively, using the `couenne` executable; AMPL also provides
   [Couenne executables](https://ampl.com/products/solvers/open-source/#couenne)),
 * GAMS (through the [GAMSlinks](https://github.com/coin-or/GAMSlinks) project),
 * [JuMP](https://github.com/jump-dev/AmplNLWriter.jl) (Julia, via AmplNLWriter),
 * [Pyomo](https://github.com/Pyomo/pyomo) (Python),
 * NEOS, in [AMPL format](https://neos-server.org/neos/solvers/minco:Couenne/AMPL.html)
   or [GAMS format](https://neos-server.org/neos/solvers/minco:Couenne/GAMS.html).

This project was initiated in 2006 within a [collaboration](http://egon.cheme.cmu.edu/ibm/page.htm)
between IBM and Carnegie Mellon University, aimed at developing algorithms for
MINLP. 

Couenne is written in C++ and is released as open source under the [Eclipse Public License 1.0](http://www.opensource.org/licenses/eclipse-1.0).

It is distributed under the auspices of the [COIN-OR Foundation](https://www.coin-or.org).

The Couenne development site is https://github.com/coin-or/Couenne.

## CITE

Paper: https://doi.org/10.1080/10556780903087124

## CURRENT BUILD STATUS

[![Windows Builds](https://github.com/coin-or/Couenne/actions/workflows/windows-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/Couenne/actions/workflows/windows-ci.yml?query=branch%3Amaster)

[![Linux and MacOS Builds](https://github.com/coin-or/Couenne/actions/workflows/linux-ci.yml/badge.svg?branch=master)](https://github.com/coin-or/Couenne/actions/workflows/linux-ci.yml?query=branch%3Amaster)

## DOWNLOAD

What follows is a quick start guide for obtaining or building
Couenne on common platforms. More detailed information is
available [here](https://coin-or.github.io/user_introduction.html).

### Docker image

There is a Docker image that provides Couenne, as well as other projects
in the [COIN-OR Optimization
Suite](https://github.com/coin-or/COIN-OR-OptimizationSuite) [here](https://hub.docker.com/repository/docker/coinor/coin-or-optimization-suite)

### Binaries

For newer releases, binaries will be made available as assets attached to
releases in Github
[here](https://github.com/coin-or/Couenne/releases). Older binaries
are archived as part of Couenne
[here](https://www.coin-or.org/download/binary/Couenne).

 * *Linux* (see https://repology.org/project/coin-or-couenne/versions for a complete listing): 
   * Fedora/Redhat/CentOS:
     ```
     $ sudo dnf install  coin-or-Couenne coin-or-Couenne-devel
     ```
   * freebsd:
     ```
     $ sudo pkg install math/couenne
     ```
 * *Windows*: The easiest way to get Couenne on Windows is to download an archive as described above.


Due to license incompatibilities, pre-compiled binaries lack some 
functionality. If binaries are not available for your platform for the latest 
version and you would like to request them to be built and posted, feel free 
to open a discussion on Github. 

### Source

Source code can be obtained either by

 * Downloading a snapshot of the source code for the latest release version of Couenne from the
 [releases](https://github.com/coin-or/Couenne/releases) page,
 * Cloning this repository from [Github](https://github.com/coin-or/Couenne), or 
 * Using the [coinbrew](https://github.com/coin-or/coinbrew) script to get the project and all dependencies (recommended, see below).   

### Dependencies

Couenne has a number of dependencies, which are detailed in
[config.yml](.coin-or/config.yml). Dependencies on other COIN-OR projects are
automatically downloaded when obtaining the source with `coinbrew`. For some
of the remaining third-party dependencies, automatic download scripts and
build wrappers are provided (and will also be automatically run for required
and recommended dependencies), while other libraries that are aeasy to obtain
must be installed using an appropriate package manager (or may come with your
OS by default). 

## BUILDING from source

These quick start instructions assume you are in a bash shell. 

### Using `coinbrew`

To download and build Couenne from source, execute the 
following on the command line. 
```
wget https://raw.githubusercontent.com/coin-or/coinbrew/master/coinbrew
chmod u+x coinbrew
./coinbrew fetch Couenne@master
./coinbrew build Couenne
```
For more detailed instructions on coinbrew, see https://coin-or.github.io/coinbrew.
The `coinbrew` script will fetch the additional projects specified in the Dependencies section of [config.yml](.coin-or/config.yml).

### Without `coinbrew` (Expert users)

 * Download the source code, e.g., by cloning the git repo https://github.com/coin-or/Couenne
 * Download and install the source code for the dependencies listed in [config.yml](.coin-or/config.yml)
 * Build the code as follows (make sure to set PKG_CONFIG_PTH to install directory for dependencies).

```
./configure -C
make
make test
make install
```

### Linear solvers for Ipopt

Couenne solves its NLP subproblems with Ipopt, which needs a sparse linear
solver. When building with `coinbrew`, [MUMPS](https://github.com/coin-or-tools/ThirdParty-Mumps)
is built by default. The HSL codes (e.g., MA27) are often faster. To use them,
obtain the source from [HSL](https://www.hsl.rl.ac.uk/ipopt/) and follow the
instructions in [ThirdParty-HSL](https://github.com/coin-or-tools/ThirdParty-HSL)
before building.

## USAGE

Couenne is run as follows:

```
couenne instance.nl
```

where `instance.nl` is an AMPL stub (`.nl`) file. Such files can be
generated from AMPL with the command "`write gfilename;`" (notice the "g"
before the file name), for example.

You may also specify a set of options to tweak the performance of Couenne.
These are read from the `couenne.opt` option file in the current directory.
A sample option file is given in [src/couenne.opt](src/couenne.opt).

## Doxygen Documentation

If you have `Doxygen` available, you can build a HTML documentation by typing

`make doxydoc` 

in the build directory. If Couenne was built via `coinbrew`, then the build
directory will be `./build/Couenne/master` by default. The doxygen documentation main file
is found at `<build-dir>/doxydoc/html/index.html`.

If you don't have `doxygen` installed locally, you can use also find the
documentation [here](http://coin-or.github.io/Couenne/Doxygen).

## DOCUMENTATION

 * [User's Manual](doc/couenne-user-manual.pdf), with explanations on most of
   the options available in Couenne
 * [List of options](doc/options.md)

## CONTRIBUTORS

 * [Pietro Belotti](http://myweb.clemson.edu/~pbelott) (Xpress Development Team, FICO, Birmingham UK)
 * [Timo Berthold](http://www.zib.de/berthold) (Xpress Development Team, FICO, Berlin)
 * [Pierre Bonami](http://pageperso.lif.univ-mrs.fr/~pierre.bonami/index.html) (IBM Cplex)
 * [Sonia Cafieri](http://www.recherche.enac.fr/~cafieri) (École Nationale de l'Aviation Civile)
 * [François Margot](http://wpweb2.tepper.cmu.edu/fmargot/index.html) (Carnegie Mellon University)
 * [Cameron Megaw](https://mthsc.clemson.edu/directory/view_person.py?person_id=218) (Clemson University)
 * [Stefan Vigerske](http://www.gams.com/~stefan/) (GAMS)
 * [Andreas Wächter](http://users.iems.northwestern.edu/~andreasw) (Northwestern University)

## ACKNOWLEDGMENTS

Credit should be given to our colleagues in the IBM/Carnegie Mellon
collaboration: Andreas, François, Pierre, Stefan, and Timo, who developed part
of Couenne, and
[Larry T. Biegler](http://www.cheme.cmu.edu/people/faculty/lb01.htm),
[Gérard Cornuéjols](http://integer.tepper.cmu.edu/),
[Ignacio E. Grossmann](http://www.cheme.cmu.edu/people/faculty/grossmann.htm), and
[Jon Lee](https://sites.google.com/site/jonleewebpage).
Each has contributed an essential part of the development of Couenne.


## Project Links

 * [Code of Conduct](https://www.coin-or.org/code-of-conduct/)
 * [COIN-OR Web Site](http://www.coin-or.org/)
 * [COIN-OR general discussion forum](https://github.com/orgs/coin-or/discussions)
 * [Couenne Discussion forum](https://github.com/coin-or/Couenne/discussions)
 * [Report a bug](https://github.com/coin-or/Couenne/issues/new)
 * [Doxygen generated documentation](http://coin-or.github.io/Couenne/Doxygen)

