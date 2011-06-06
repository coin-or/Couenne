OPT = -O3
PROF =

ifeq ($(HOSTNAME),octa.tepper.cmu.edu)
OPT = -O3
PROF =
endif
ifeq ($(HOSTNAME),pbonami-lx.tepper.cmu.edu)
OPT = -O3
PROF =
endif

LIBS = 

SOLVER = CPLEX
# Possibilites: CPLEX or CLP; exactly one must be set

DEFS := -D$(SOLVER) \
	-DF_RES -DSHORT_F_RES \
	-DaaTRACE -DaaTRACE_EIGENVALUES -DaaTRACE_USED_EIGENVECTORS -DaaTRACE_ALL -DaaTRACE_RLT_MPS \
	-DaaTRACE_CUT_TIME -DaaTRACE_SPARSIFY_CUT_DEPTH \
	-DaaTRACE_HEUR \
	-DaaDETAILED_SPARSIFY_COMPARISON \
	-DaaCHECK_QUADRATIC_CUTS \
	-DaaTRACE_DISJUNCTIVE_CUTS \
	-DaaTRACE_ORTHOCUT \
	\
	-DaaUSE_INTERIOR_POINT \
	\
	-DRLT_CUTS -DRLT_SEPARATION \
	-DINCL_ORIG -DaaSCALE_EIGENV -DONLY_NEG_EIGENV -DaaONLY_MOST_NEG \
	-DaaNORMALIZE_SPARSE_CUTS -DaaSPARSIFY2 -DaaSPARSIFY \
	-DSPARSIFY_OLD_DELTA=0.50        -DSPARSIFY_NEW_DELTA=0.50 \
	-DSPARSYFY_OLD_NZ_THRESHOLD=0.50 -DSPARSYFY_NEW_NZ_THRESHOLD=0.70 \
	-DWISE_SPARSIFY -DWISE_SPARSIFY_GAP=0.0001 \
	-DSPARSIFY_MINOR_SDP_CUTS \
	-DaaSPARSIFY_REMOVE_DUPLICATES -DSPARSIFY_MAX_CARD=10000\
	-DaaZERO_LB_DIAGONAL \
	-DaaSLIDING \
	-DaaORTHOCUT \
	-DaaDISJUNCTIVE_CUTS \
	-DaaLINQUAD_CUTS \
	-DaaLINQUAD_BOUNDS_CUTS \
	\
	-DaaADD_ONLY_PERC_DEEPEST_SPARSIFY_CUTS=0.20 \
	-DaaADD_RAND_PERC_SPARSIFY_CUTS=0.20 \
	-DMAX_CUTS_PER_ITER=1000 \
	-DaaRAND_CUT_ADD=0.33 \
	-DDELETE_INACTIVE_CUTS=1e-5 -DaaFANCY_CUT_ACTIVITY=3 \
	-DCUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC=0.0001 -DCUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER=3\
	\
	-DHEUR_XXT -DHEUR_MIN_MATRIX_NORM -DHEUR_IMPROVE_SOLUTION\
	\
	-DEXIT_ON_ITER=1000 -DEXIT_ON_TIME=300.0 \
	-DaaEXIT_IMPROVEMENT_ITER=10 -DaaEXIT_TAILING_OFF_VALUE=0.001 -DaaEXIT_TAILING_OFF_PERC=0.01 \
	\
	-DCHECK -DaaWRITE_ITER_LP

# Possibilities:
#
# Trace:
#                -DTRACE: trace
#                -DTRACE_EV: trace eigenvalues of the current solution
#                -DTRACE_ALL: detailed trace
#                -DTRACE_RLT_MPS: write the initial RLT formulation in rlt.mps
#                -DTRACE_HEUR: trace heuristics
#                -DF_RES: Write a detailed solver log in f_res.xxx
#                -DSHORT_F_RES: Append a solution summary line in short_f_res.xxx
#
# Solver options:
#                -USE_INTERIOR_POINT: uses CPLEX interior point solver baropt with disabled crossover
#                                     for solving the initial relaxation
#
# Alg. options:
#		 -DRLT_CUTS: add RLT cuts to original formulation
#		 -DRLT_SEPARATION: does not include ALL the RLT alltogether but separates them (RLT_CUTS must be enabled)
#		 -DZERO_LB_DIAGONAL: set the lower bound of X_ii to 0
#                -DINCL_ORIG: include cut from original eigenvector
#                -DSCALE_EIGENV: Multiply eigenvector by sqrt(n) before
#                                sparsifying
#                -DONLY_NEG_EIGENV: use only eigenvectors for negative 
#                                   eigenvalues
#                -DONLY_MOST_NEG: use only eigenvector for most negative 
#                                 eigenvalue
#                -DSPARSIFY: zero some components of eigenvectors before
#                            generating cut
#                -DSPARSIFY_MINOR_SDP_CUTS: separate SDP cuts from minor submatrices 
#                                             of {{1,x},{x,X}} induced by sparsified cuts
#                -DWISE_SPARSIFY: recomputes eigendecomposition each time we zero some component
#                -DSLIDING: generate eigencuts, x := x + \lambda M and repeat
#
# Cut management:
#                -DMAX_CUTS_PER_ITER=<value>: adds at most <value> new cuts per iteration
#                -DRAND_CUT_ADD=<value>: add generated cuts with probability <value>
#                -DADD_ONLY_PERC_DEEPEST_CUTS=<value>: add only the (100*<value>)% of the cuts 
#                                                      (ordered by depth)
#                -DDELETE_INACTIVE_CUTS=<value>: remove all cuts whose dual vars is 
#                                                smaller than <value>
#                -DFANCY_CUT_ACTIVITY=<value>: removes a cut only only if it has been inactive 
#                                              in the last <value> iterations
#		 -DCUT_DELETION_OBJ_PERC_IMPROVEMENT_PERC<value>
#		    -DCUT_DELETION_OBJ_PERC_IMPROVEMENT_ITER<iter>
#		                               removes cuts only if the objective 
#					       function increment was MORE than (100*<value>)% 
#					       from the past <iter> iterations
#
# Heuristics:
#                -DHEUR_XXT: considers the solution (x,xx^t) where x is the current
#                            outer-approximation LP solution (note: it might be unfeasible for QCQP)
#                -DHEUR_MIN_MATRIX_NORM: computes a solution "close" to the outer-approximation
#                                        solution (x,X)
#                -DHEUR_IMPROVE_SOLUTION: if the problem has y_i variables it might be possible
#                                         to improve the current heuristic solution by fixing
#                                         the variables x_i and X_ij and optimizing over y_i
# Exit:
#                -DEXIT_ON_ITER=<value>: exit after <value> iterations
#                -DEXIT_ON_TIME=<value>: exit after <value> seconds
#                -DEXIT_IMPROVEMENT_ITER=<n>
#                  -DEXIT_TAILING_OFF_VALUE=<value>: exits if the solution value improved 
#                                                    LESS than <value> in the last <n> 
#                                                    iterations
#                  -DEXIT_TAILING_OFF_VALUE=<value>: exits if the solution value improved 
#                                                    LESS than (100*<value>)% from the previous 
#                                                    <n>th iteration
#
# Debug:
#                -DCHECK: misc. runtime checks
#                -DWRITE_ITER_LP: write lastiter.lp at the end of each iteration
#

HOMEBIN=.

CXX := g++
CC  := gcc

EXTRAS =
OSICPXLIB = -lOsiCpx -lcplex 


ifeq ($(HOSTNAME),pbonami-lx.tepper.cmu.edu)
#  USRINC=${HOME}/.usr/include/
#  USRLIB=$(HOME)/.usr/lib/
  USRINC = /home/aqualizz/AsxCOIN/include/
  USRLIB = /home/aqualizz/AsxCOIN/lib/
  EXTRAS += /usr/lib/liblapack.a /usr/lib/libblas.a 
  CPLEXINC= -I/usr/ilog/cplex110/include/ilcplex
  CPLEXLIB= -L/usr/ilog/cplex110/lib/x86_rhel4.0_3.4/static_pic/
  LIBS += -lAsxBmqs -lAsxBp -lAsxCpx -lAsxUtils
endif
ifeq ($(HOME),/home/fmargot)
  USRINC = ${ENV_COIN_BUILD}/include/
  USRLIB = $(ENV_COIN_BUILD)/lib/
  EXTRAS += /usr/lib/liblapack.a /usr/lib/libblas.a
  CPLEXINC = -I$(ENV_CPLEX_H)
  CPLEXLIB = -L$(ENV_CPLEX_LIB)
endif

ifeq ($(HOME),/home/andrea)
#  USRINC = ${ENV_COIN_BUILD}/include/
#  USRLIB = ${ENV_COIN_BUILD}/lib/
  USRINC = /home/andrea/AsxCOIN/include/
  USRLIB = /home/andrea/AsxCOIN/lib/
  EXTRAS += /usr/lib/liblapack.a 
#  CPLEXINC = -I$(ENV_CPLEX_H)
#  CPLEXLIB = -L$(ENV_CPLEX_LIB)
  CPLEXINC = -I/usr/ilog/cplex111/include/ilcplex/
  CPLEXLIB = -L/usr/ilog/cplex111/lib/x86_debian4.0_4.1/static_pic/
#  CPLEXINC = -I/usr/ilog/cplex110/include/ilcplex/
#  CPLEXLIB = -L/usr/ilog/cplex110/lib/x86_rhel4.0_3.4/static_pic/
  LIBS += -lAsxBmqs -lAsxBp -lAsxCpx -lAsxUtils
else ifeq ($(HOSTNAME),octa.tepper.cmu.edu)
#  USRINC=${HOME}/.usr/include/
#  USRLIB=$(HOME)/.usr/lib/
  USRINC = /home/aqualizz/AsxCOIN/include/
  USRLIB = /home/aqualizz/AsxCOIN/lib/
  EXTRAS += /usr/lib64/liblapack.so /usr/lib64/libblas.so 
  CPLEXINC= -I/usr/ilog/cplex111/include/ilcplex
  CPLEXLIB= -L/usr/ilog/cplex111/lib/x86-64_debian4.0_4.1/static_pic/
  LIBS += -lAsxBmqs -lAsxBp -lAsxCpx -lAsxUtils
endif

# include paths
CGLINC = -I${USRINC}
OSIINC = -I${USRINC}
COINUTILSINC = -I${USRINC}

# library paths
CGLLIB = ${USRLIB}
OSILIB = ${USRLIB}
COINUTILSLIB = ${USRLIB}

#



CFLAGS = -Wall $(OPT) $(PROF) $(DEFS) \
	${CGLINC} ${COINUTILSINC} ${OSIINC} ${CPLEXINC} -I.

LFLAGS = -Wall $(OPT) $(PROF) \
	-L${CGLLIB} -L${COINUTILSLIB} -L${OSILIB} ${CPLEXLIB} \

ifeq ($(HOSTNAME),farkas.tepper.cmu.edu)
  LFLAGS += -L${ENV_G2C_LIB} 
endif

OBJ =	populate.o misc_util.o Heuristics.o CutGen.o dsyevx_wrapper.o tracer.o quadratic_cuts_check.o linquad_cuts.o orthocut.o disjunctive_cuts.o rlt_cuts.o

#LIBS = -lstdc++ -lm -lgfortran \
#       -lCgl -lOsi -lOsiClp $(OSICPXLIB) \
#       -lpthread -lClp -lCoinUtils -lg2c
LIBS += -lstdc++ -lm -lgfortran \
       -lCgl -lOsi -lOsiClp $(OSICPXLIB) \
       -lpthread -lClp -lCoinUtils -llapack -lblas

HEADER= sdpcuts.hpp populate.hpp misc_util.hpp Heuristics.hpp CutGen.hpp dsyevx_wrapper.hpp tracer.hpp quadratic_cuts_check.hpp linquad_cuts.hpp orthocut.hpp disjunctive_cuts.hpp rlt_cuts.hpp

all: ${HOMEBIN}/sdpcuts ${HOMEBIN}/sedumi ${HOMEBIN}/report ${HOMEBIN}/report2 ${HOMEBIN}/tracerreport

${HOMEBIN}/sdpcuts: sdpcuts.o ${OBJ} ${HEADER}
	$(CC) $(LFLAGS) \
		sdpcuts.o ${OBJ} \
		${EXTRAS} \
		${LIBS} \
		-o ${HOMEBIN}/sdpcuts

${HOMEBIN}/sedumi: sedumi.o populate.o misc_util.o ${HEADER}
	$(CC) $(LFLAGS) \
		sedumi.o  populate.o misc_util.o \
		${EXTRAS} \
		${LIBS} \
		-o ${HOMEBIN}/sedumi

${HOMEBIN}/report: report.o
	$(CC) $(LFLAGS) \
		-lstdc++ report.o \
		-o ${HOMEBIN}/report

${HOMEBIN}/report2: report2.o misc_util.o
	$(CC) $(LFLAGS) \
		${LIBS} \
		-lstdc++ report2.o misc_util.o\
		-o ${HOMEBIN}/report2


${HOMEBIN}/tracerreport: misc_util.o tracerreport.o
	$(CC) $(LFLAGS) \
		${LIBS} \
		-lstdc++ misc_util.o tracerreport.o \
		-o ${HOMEBIN}/tracerreport

%.o: %.cpp Makefile
	@echo $(HOSTNAME)

	$(CXX) $(CFLAGS) -c $< 

clean:
	rm -f $(OBJ) sdpcuts.o sdpcuts report.o report

cleano:
	rm -f *.o

clean~:
	rm -f *.c~ *.h~ *.cpp~ *.hpp~ Makefile~

cleanlogs:
	rm -f f_res.xxx short_f_res.xxx detailedreport.txt globalreport.txt

