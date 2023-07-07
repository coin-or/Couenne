# Couenne

Couenne (Convex Over and Under ENvelopes for Nonlinear Estimation) is a branch&bound algorithm to solve Mixed-Integer Nonlinear Programming (MINLP) problems of the form:

&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; min <i>f</i><sub>0</sub>(<i>x</i>,<i>y</i>)<br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     <i>f</i><sub>i</sub>(<i>x</i>,<i>y</i>) &le; 0&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;<i>i</i>=1,2..., <i>m</i><br/>
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;     <i>x</i> in R<sup>n</sup>,
                                   <i>y</i> in Z<sup>p</sup>

where all _f_<sub>i</sub>(_x_,_y_) are, in general, nonlinear functions. 

Couenne aims at finding global optima of nonconvex MINLPs. It implements linearization, bound reduction, and branching methods within a branch-and-bound framework. Its main components are:

 * an expression library;
 * separation of linearization cuts;
 * branching rules;
 * bound tightening methods.

It is distributed on [COIN-OR](http://www.coin-or.org) under the [Eclipse Public License](http://www.opensource.org/licenses/eclipse-1.0) (EPL). The EPL is a license approved by the Open Source Initiative (OSI), thus Couenne is OSI Certified Open Source Software.


## Download, installation and usage

Couenne is found on the COIN-OR [project page](https://github.com/coin-or/Couenne). It can be downloaded with [Subversion](http://subversion.tigris.org) -- see also some [instructions](https://projects.coin-or.org/BuildTools/wiki/user-download#ObtainingtheCodeUsingSubversion) on using svn. Run the command:

```
svn co https://projects.coin-or.org/svn/Couenne/stable/0.5 Couenne
```

to get the source code of the stable version. Before building and installing Couenne, some third party packages are needed. These cannot be downloaded from COIN-OR, and have to be obtained independently.
These packages are: ASL, Blas, Lapack, and HSL or MUMPS.
The user is referred to the `INSTALL` file in each of the subdirectories of Couenne/ThirdParty for instructions on how to obtain them.

In general, the stable version of Couenne is subject to slight changes such as bug fixes. In order to be up-to-date with such changes, you may run the command `svn update` within the `Couenne` directory.
All releases are also available as archive at https://www.coin-or.org/download/source/Couenne/. A release is not subject to change.

To install Couenne, we refer to general installation [instructions](https://projects.coin-or.org/BuildTools/) for COIN-OR projects.
We also suggest the excellent [Ipopt compilation hints](https://projects.coin-or.org/Ipopt/wiki/CompilationHints) page for compiling on non-Linux systems, such as Mac and Windows.
The impatient may want to issue the following commands:

```
cd Couenne
cd ThirdParty      # Read INSTALL.* file in each subdirectory and get third party software
cd ..
mkdir build
cd build
../configure -C
make
make install
```

The above commands place Couenne in the `Couenne/build/bin/` directory, libraries in `Couenne/build/lib/`, and include files in `Couenne/build/include/`. An alternative directory can be specified with
the `--prefix` option of configure. For instance, when replacing "`../configure -C`" above with

```
../configure -C --prefix=/usr/local
```

the Couenne executable will be installed in `/usr/local/bin/`, the libraries in `/usr/local/lib/`, and the include files in `/usr/local/include/`.
Couenne is run as follows:

```
couenne instance.nl
```

where `instance.nl` is an AMPL stub (`.nl`) file. Such files can be generated from AMPL with the command "`write gfilename;`" (notice the "g" before the file name), for example.

You may also specify a set of options to tweak the performance of Couenne. These are found in the `couenne.opt` option file. A sample option file is given in the Couenne/src/ directory.


## Documentation

A user [manual](doc/couenne-user-manual.pdf) is available, with explanations on most options available in Couenne. Doxygen documentation is also available, and it can be generated by running

```
make doxydoc
```

from the same `build/` directory where you ran configure, make, and make
install. Documentation in both html and LaTeX format can be found in
the Doc/ subdirectory. Fire up your browser and take a look at
Doc/html/index.html for documentation of Couenne.



## Resources and links

Couenne is maintained by [Pietro Belotti](https://belotti.faculty.polimi.it).

Web page: [https://www.github.com/coin-or/Couenne](https://www.github.com/coin-or/Couenne)

Dependencies: 
[CoinUtils](https://github.com/coin-or/CoinUtils), 
[Cbc](https://github.com/coin-or/Cbc), 
[Cgl](https://github.com/coin-or/Cgl), 
[Clp](https://github.com/coin-or/Clp), 
[Ipopt](https://github.com/coin-or/Ipopt), and 
[Osi](https://github.com/coin-or/Osi) (from COIN-OR);
ASL (the [Ampl](http://www.ampl.com) Solver Library),
Lapack,
Blas,
HSL,
MUMPS,
[SCIP](http://scip.zib.de), and
[SoPlex](http://soplex.zib.de).

External resources:  [COIN-OR](http://www.coin-or.org), 
[Eclipse Public License](http://www.opensource.org/licenses/eclipse-1.0).



## Report a bug, contribute to Couenne

As an open-source code, contributions to Couenne are welcome. To submit a contribution to Couenne, please follow the [COIN-OR guidelines](http://www.coin-or.org/contributions.html).

In order to report a bug, use the [issue system](https://github.com/coin-or/Couenne/issues).

In order to ensure that your issue is addressed in a timely fashion, please try to be as exhaustive as you can in the bug report, for instance by reporting what version of Couenne you have downloaded and what operating system you are using, and again by attaching the model/data files with which the crash occurred.


## Contributors

 * [Pietro Belotti](http://myweb.clemson.edu/~pbelott) (Xpress Development Team, FICO, Birmingham UK)
 * [Timo Berthold](http://www.zib.de/berthold) (Xpress Development Team, FICO, Berlin)
 * [Pierre Bonami](http://pageperso.lif.univ-mrs.fr/~pierre.bonami/index.html) (IBM Cplex)
 * [Sonia Cafieri](http://www.recherche.enac.fr/~cafieri) (École Nationale de l'Aviation Civile)
 * [François Margot](http://wpweb2.tepper.cmu.edu/fmargot/index.html) (Carnegie Mellon University)
 * [Cameron Megaw](https://mthsc.clemson.edu/directory/view_person.py?person_id=218) (Clemson University)
 * [Stefan Vigerske](http://www.gams.com/~stefan/) (GAMS)
 * [Andreas Wächter](http://users.iems.northwestern.edu/~andreasw) (Northwestern University)



## Acknowledgments

This project was initiated in 2006 within a [collaboration](http://egon.cheme.cmu.edu/ibm/page.htm) between IBM and Carnegie Mellon University, aimed at developing algorithms for MINLP.

Credit should be given to our colleagues in this collaboration: Andreas, François, Pierre, Stefan, and Timo, who developed part of Couenne, and
[Larry T. Biegler](http://www.cheme.cmu.edu/people/faculty/lb01.htm), 
[Gérard Cornuéjols](http://integer.tepper.cmu.edu/),
[Ignacio E. Grossmann](http://www.cheme.cmu.edu/people/faculty/grossmann.htm), and
[Jon Lee](https://sites.google.com/site/jonleewebpage).
Each has contributed an essential part of the development of Couenne.


## Project Links

 * [Couenne executables](http://ampl.com/products/solvers/open-source/#couenne) provided by AMPL
 * Couenne is available through the [JuMP](https://github.com/JuliaOpt/JuMP.jl) modeling language
 * [COIN-OR Initiative](http://www.coin-or.org/)
 * [IBM CMU Open Source MINLP project](http://egon.cheme.cmu.edu/ibm/page.htm)  
 * [Report a bug](https://github.com/coin-or/Couenne/issues)
 * [NEOS server for Optimization](https://neos-server.org/neos/): You can submit a MINLP problem for solution with Couenne in [AMPL format](https://neos-server.org/neos/solvers/minco:Couenne/AMPL.html) or [GAMS format](https://neos-server.org/neos/solvers/minco:Couenne/GAMS.html)


--------

## Options for BonCouenne


### Linearization options

```
convexification_cuts <num>
```

Specify the frequency (in terms of nodes) at which linearization cuts are generated. Default: 1.
If 0, linearization cuts are never separated.

```
convexification_points <num>
```

Specify the number of points at which to convexify. Default: 1.

```
violated_cuts_only <yes|no>
```

If set to yes (default), only violated convexification cuts will be added.


```
art_lower <num>
```

Set artificial lower bound (for minimization problems),
useful when a lower bound is known or for testing purposes.
Default value is -10<sup>50</sup>.

```
opt_window <num>
```

Multiplier for restricting variable bounds around known optimum (to be read from file with 
method CouenneProblem::readOptimum()). If the optimal value x,,i,, of the i-th variable is known,
before starting Couenne its bounds will be intersected with interval [x<sub>i</sub>-K(1+|x<sub>i</sub>|),x<sub>i</sub>+K(1+|x<sub>i</sub>|)],
where K is the value of the option. Default value is infinity.


```
use_quadratic <yes|no>
```

Use quadratic expressions and related exprQuad class. Still in testing, so default is "no".


### Bound tightening options

```
feasibility_bt <yes|no>
```

Use feasibility-based bound tightening (strongly recommended).
Default value is "yes".


```
optimality_bt
```

Optimality-based (expensive) bound tightening. Only recommended for problems with few variables 
and/or at the initial nodes of the B&B tree. Default is "no". If set to "yes", we recommend 
to couple it with a value of log_num_obbt_per_level of 0 (see below).

```
log_num_obbt_per_level <num>
```

Specify the frequency (in terms of nodes) for optimality-based bound tightening. Default is 0.

 * If -1, apply at every node (expensive!).
 * If 0, apply at root node only.
 * If k>0, apply with probability 2<sup>(k - level)</sup>, level being the current depth of the B&B tree.

```
aggressive_fbbt <yes|no"
```

Aggressive feasibility-based bound tightening (to use with NLP points). Default value is "yes". 
This is also computationally expensive.

```
log_num_abt_per_level <num>
```

Specify the frequency (in terms of nodes) for aggressive bound tightening (similar to log_num_obbt_per_level).

 * If -1, apply at every node (expensive!);
 * If 0, apply at root node only;
 * If k>0, apply with probability 2<sup>(k - level)</sup>, level being the current depth of the B&B tree.



### Branching options

```
branch_fbbt <yes|no>
```

Apply bound tightening before branching. default: yes

```
branch_conv_cuts <yes|no>
```

Apply convexification cuts before branching (not active yet). Default: no.

```
branch_pt_select <string>
```

Chooses branching point selection strategy. Possible values are

 * "lp-clamped": LP point clamped in [k,1-k] of the bound intervals (k defined by lp_clamp);
 * "lp-central": LP point if within [k,1-k] of the bound intervals, middle point otherwise (k defined by  branch_lp_clamp);
 * "balanced":   minimizes max distance from curve to convexification;
 * "min-area":   minimizes total area of the two convexifications;
 * "mid-point":  convex combination of current point and mid point;
 * "no-branch":  do not branch, return null infeasibility; for testing purposes only.

Default is mid-point.

```
branch_lp_clamp <num>
```

Defines a threshold for selecting an LP point as the branching point;
<num> is between 0 and 0.5 and defaults to 0.2. Suppose variable x,,i,, with bounds [l<sub>i</sub>,u<sub>i</sub>]
is chosen for branching. If the lp-central or lp-clamp strategies 
are selected, the branching point is projected into the interval [L<sub>i</sub>,U<sub>i</sub>] with
L<sub>i</sub> = l<sub>i</sub>+ a(u<sub>i</sub> - l<sub>i</sub>) and U<sub>i</sub> = l<sub>i</sub>+ (1-a)(u<sub>i</sub> - l<sub>i</sub>).

```
branch_midpoint_alpha <num>
```

Defines convex combination of mid point and current LP point: branching point will be
alpha x<sub>i</sub> + (1-alpha) (l<sub>i</sub>+u<sub>i</sub>)/2. Default value is 0.25.

Options `branch_pt_select` and `branch_lp_clamp` above are
also available for the following set of operators, and are applied to 
each of these operators independently:
"prod", "div", "exp", "log", "trig", "pow",  "negpow", "sqr", "cube".

For instance, the following settings: 

```
branch_pt_select balanced

branch_pt_select_prod lp-clamp
branch_lp-clamp_prod 0.15

branch_pt_select_log lp-central
branch_lp-clamp_log 0.1
```

specify balanced strategy for all operators except 
products and logarithms, lp-clamp with parameter 0.15 for products and
lp-central with parameter 0.1 for logarithms.



### Upper bounding options

```
local_optimization_heuristic <yes|no>
```

Search for local solutions of NLPs. Default: yes.

```
log_num_local_optimization_per_level <num>
```

Specify the logarithm of the number of local optimizations to perform
on average for each level of given depth of the tree.

If equal to -1, solve as many nlp's at the nodes for each level of the tree.

Nodes are randomly selected. If for a given level there are less 
nodes than this number nlp, are solved for every nodes.
For example, if parameter is 8, nlp's are solved for all node until level 8,
then for half the node at level 9, 1/4 at level 10.


```
art_cutoff <num>
```

Set artificial cutoff useful when a feasible solution is known or for testing purposes.
Default value is 10<sup>50</sup>.


```
feas_tolerance
```

This is a feasibility tolerance for candidate feasible solutions. 
Default value is 10<sup>-7</sup>.


<!--
### TODO

```
pseudocost_mult
```


```
pseudocost_mult_lp
```


```
enable_sos
```


```
branch_fbbt
```


```
branch_conv_cuts
```


```
branch_pt_select
```


```
branch_midpoint_alpha
```


```
branch_lp_clamp
```


```
cont_var_priority
```


```
red_cost_branching
```


```
convexification_cuts
```


```
check_lp
```


```
local_optimization_heuristic
```


```
log_num_local_optimization_per_level
```


```
convexification_type
```


```
convexification_points
```


```
violated_cuts_only
```


```
art_cutoff
```


```
opt_window
```


```
feas_tolerance
```


```
feasibility_bt
```


```
use_quadratic
```


```
optimality_bt
```


```
log_num_obbt_per_level
```


```
aggressive_fbbt
```


```
log_num_abt_per_level
```


```
art_lower
```


```
branching_object
```
-->
