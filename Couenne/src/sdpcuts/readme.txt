

CutGen.cpp/hpp			SDP cut generators + sparse SDP cut generators
disjunctive_cuts.cpp/hpp	Interface with Anureet Disjunctive Cut generator. We never really used it and the code needs to be reviewed. Not used in our experiments [IGNORE]
dsyevx_wrapper.cpp/hpp		wrappers for the Lapack eigendecomposition calls
Heuristics.cpp/hpp		two heuristics to obtain solutions in the original space (x,y), given one in the higher dimensional space (x,X,y) - not used in our experiments [IGNORE]
linquad_cuts.cpp/hpp		generates linear cuts that approximate the quadratic constraint X_ii >= x_i^2; experiments showed that these cuts are not practically useful [IGNORE]
quadratic_cuts_check.cpp/hpp	given a higher dimensional solution (x,X,y) checks if quadratic constraints X_ii >= x_i^2 are violated. [IGNORE]
misc_util.cpp/hpp		utility functions
orthocut.cpp/hpp		not really useful - it is shown to be dominated by standard sdp cuts [IGNORE]
OsiXxxSolverInterface.hpp	extension to OsiCpxSolverInterface
populate.cpp/hpp		reads a problem in mps format and initializes the working structures
report(2).cpp/hpp		used to generate a report when comparing two executions with different parameters given the f_res.xxx output files [IGNORE]
rltcuts.cpp/hpp			generates RLT inequalities
sdpcuts.cpp/hpp			main algorithm. Where the "main" is
sdpsol.cpp/hpp/sedumi.cpp	given an MPS problem reads it and produces a .mat (Matlab) file ready to be used with SeDuMi or CSDP [IGNORE]
tracer.cpp/hpp			class used to collect info from all components of the algorithm. Generates a detailed report and a global report (with less info)
tracerreport.cpp 		separate program. Compares two globalreports from two different executions of the algorithm.



