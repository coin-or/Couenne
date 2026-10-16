# Options for BonCouenne


## Linearization options

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


## Bound tightening options

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



## Branching options

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



## Upper bounding options

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
## TODO

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
