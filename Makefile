##############################################################################
##############################################################################
#                                                                            #
# This file is part of the SYMPHONY Branch, Cut, and Price Library.          #
#                                                                            #
# SYMPHONY was jointly developed by Ted Ralphs (tkralphs@lehigh.edu) and     #
# Laci Ladanyi (ladanyi@us.ibm.com).                                         #
#                                                                            #
# (c) Copyright 2000, 2001, 2002  Ted Ralphs. All Rights Reserved.           #
#                                                                            #
# This software is licensed under the Common Public License. Please see      #
# accompanying file for terms.                                               #
#                                                                            #
##############################################################################
##############################################################################

##############################################################################
##############################################################################
# These entries define the various architectures that are supported as
# determined by the variable ARCH. If your arcitecture is not listed, simply
# add it below and then search for the places where there are
# architecture-specific variables set and make sure to set them properly for
# your specific architecture. For each architecture, there will be three
# subdirectories created called $(USERROOT)/bin.$(ARCH),
# $(USERROOT)/dep.$(ARCH) and ${USERROOT)/objects.$(ARCH) where the
# corresponding objects, binaries, and dependencies for each architecture type
# will reside.
##############################################################################
##############################################################################

##############################################################################
# !!!!!!!!! User must set this variable to the proper architecture !!!!!!!!!!!
#
# CURRENT OPTIONS: LINUX, RS6K, SUN4SOL2, SUNMP, X86SOL2, ALPHA
#
# Note that these architecture variable names are the ones used to define
# the various architectures recognized by the Parallel Virtual Machine (PVM)
# library. 
##############################################################################

ARCH = LINUX

##############################################################################
# If you have PVM installed, this will set the variable ARCH automatically.
##############################################################################

#ARCH=$(PVM_ARCH)

MAKE = make

AWK = awk

AR = ar -r

RANLIB = ranlib

##############################################################################
# The ROOT environement variable specifies the root directory for the source
# code. If this file is not in the SYMPHONY root directory, change this
# variable to the correct path.
##############################################################################

ROOT = .

##############################################################################
# USERROOT is the name of the user's source directory. $(ROOT)/VRP is the
# proper setting for the sample application, a VRP and TSP solver.
##############################################################################

USERROOT = $(ROOT)/Vrp

##############################################################################
##############################################################################
# LP solver dependent definitions
##############################################################################
##############################################################################

##############################################################################
##############################################################################
#You must define an LP solver in order to use the software. However, by
#default, this option is set to "NONE," which results in an error.
##############################################################################
##############################################################################

LP_SOLVER = NONE

##############################################################################
# If you want something other than CPLEX or OSL, add it here.
##############################################################################

##############################################################################
# OSL definitions
##############################################################################

#Uncomment the line below if you want to use OSL.
LP_SOLVER = OSL

#Set the paths and the name of the library
ifeq ($(LP_SOLVER),OSL)
       LPINCDIR = -I/usr/local/include
       LPLDFLAGS = -L/usr/local/lib
       LPLIB = -losl
endif

##############################################################################
# CPLEX definitions
##############################################################################

#Uncomment the line below if you want to use CPLEX.
#LP_SOLVER = CPLEX

ifeq ($(LP_SOLVER),CPLEX)
	LPINCDIR = -I/usr/local/include
	LPLDFLAGS = -L/usr/local/lib
	LPLIB = -lcplex
endif

##############################################################################
##############################################################################
# COMM Protocol definitions
##############################################################################
##############################################################################

##############################################################################
# If you want something other than PVM, add it here. 
##############################################################################

##############################################################################
# Setting COMM_PROTOCOL to NONE will allow compilation without linking the 
# PVM library, but only works if a single executable is produced (either 
# sequential or shared-memory). Change the COMM_PROTOCOL if you want to
# compile a distributed version of the code.
##############################################################################

COMM_PROTOCOL = NONE
#COMM_PROTOCOL = PVM

#Set the paths for PVM
ifeq ($(COMM_PROTOCOL),PVM)
	COMMINCDIR = -I$(PVM_ROOT)/include
	COMMLDFLAGS = -L$(PVM_ROOT)/lib/$(PVM_ARCH)
	COMMLIBS = -lgpvm3 -lpvm3
endif

##############################################################################
##############################################################################
# Set the compiler -- If an OpenMP compliant copiler is used, then a shared
#	memory version of the code will result as long as at least 
#	COMPILE_IN_LP is set to TRUE (see comments above).
##############################################################################
##############################################################################

CC = gcc

##############################################################################
# Set the optimization level
# NOTE: if it is set to "-g" then everything is compiled debugging enabled.
#       if it is set to "-O" then it'll be reset to the highest possible opt
#       otherwise it'll left alone 
##############################################################################

OPT = -O

##############################################################################
##############################################################################
# These options are for configuring the modules and have the following
# meanings:
#
# COMPILE_IN_CG: If set to true, then the CG function will be called directly 
#	from each LP solver instead of running as a separate executable. Note 
#	that the parameter "use_cg" should be set to FALSE (the default) if 
#	this option is set. The executable containing the LP solver will have 
#	the suffix _cg added to it to denote the inclusion of the cut generator
# 	function.
# COMPILE_IN_CP: As above, if this flag is set, then the cut pool resides in
#	the LP solver and the pool is scanned directly from there. Note that
#	if this option is chosen when multiple LP processes are running, then
#	they will all have their own cut pool. The executable containing the 
#	LP solver will have the suffix _cp added to it to denote the inclusion
#	of the cut generator function.
# COMPILE_IN_LP: If this flag is set, the LP solver will be called directly 
#	from the tree manager. Note that this necessarily implies that there
#	only be one LP solver. This DOES NOT automatically imply that
#	the cut generator and/or cut pool will be compiled in. The tree 
#	manager executable name will have the appropriate suffix added to it
#	to denote the inclusion of the LP solver function.
# COMPILE_IN_TM: If this flag is set, the tree manager function will be 
#	compiled directly from the master module instead of running as a 
#	separate executable. This DOES NOT imply that the LP, cut generator 
#	or cut pool functions will be compiled in. The master executable
#	name will contain a suffix indicating what functions are compiled in.
#
# Note that if an OpenMP compliant compiler is used to make a shared-memory
#	version of the code, then it is recommended to set all these
#	variables to TRUE. However, it should work with only COMPILE_IN_LP
#	set to TRUE.
##############################################################################
##############################################################################

COMPILE_IN_CG = TRUE
COMPILE_IN_CP = TRUE
COMPILE_IN_LP = TRUE
COMPILE_IN_TM = TRUE

##############################################################################
##############################################################################
# A bunch of other compile-time options
##############################################################################
##############################################################################

#__BEGIN_EXPERIMENTAL_SECTION__#
##############################################################################
# DECOMP related stuff 
##############################################################################

DECOMP = FALSE

#___END_EXPERIMENTAL_SECTION___#
##############################################################################
# Option to only process the root node (for testing root lower bounds)
##############################################################################

ROOT_NODE_ONLY = FALSE

##############################################################################
# ccmalloc related flags
##############################################################################

CCMALLOC = ccmalloc

######################################################################
# Whether to compile in the fractional branching option
######################################################################

COMPILE_FRAC_BRANCHING = FALSE

#######################################################################
# Whether to perform additional sanity checks (for debugging purposes)
#######################################################################

DO_TESTS = FALSE

#######################################################################
# More testing ....
#######################################################################

TM_BASIS_TESTS = FALSE

#######################################################################
# Additional debugging options 
#######################################################################

TRACE_PATH = FALSE
CHECK_CUT_VALIDITY = FALSE
CHECK_LP = FALSE

#######################################################################
# Additional statistics
#######################################################################

STATISTICS = FALSE

##############################################################################
# Some experimental pseudo-cost branching stuff
##############################################################################

PSEUDO_COSTS = FALSE

##############################################################################
# Path to the X11 include directory and library files
##############################################################################

X11INCDIR  =
X11LDFLAGS =

##############################################################################
# Other variables
##############################################################################

MACH_DEP  =
SYSLIBS   = 
GCCLIBDIR =

##############################################################################
##############################################################################
# OS dependent flags, paths, libraries
# Set separate variable values for each architecture here
##############################################################################
##############################################################################

##############################################################################
# LINUX Definitions
##############################################################################

ifeq ($(ARCH),LINUX)
	X11LDFLAGS = -L/usr/X11R6/lib
	ifeq ($(LP_SOLVER),CPLEX)
	   LPSOLVER_DEFS = -DSYSFREEUNIX
	endif
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM
	SYSLIBS = -lpthread #-lefence
endif

##############################################################################
# RS6K Definitions
##############################################################################

ifeq ($(ARCH),RS6K)
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM
	SYSLIBS = -lbsd
	ifeq ($(ARCH),RS6KMP)
	   SYSLIBS += -lpthreads
	endif
endif

##############################################################################
# Sun Sparc Solaris Definitions
##############################################################################

ifeq ($(ARCH),SUN4SOL2)
	X11LDFLAGS = -L/usr/local/X11/lib -R/usr/local/X11/lib
	ifeq ($(LP_SOLVER),CPLEX)
	   LPSOLVER_DEFS = -DSYSGNUSOLARIS
	endif
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM 
	SYSLIBS = -lsocket -lnsl
endif

##############################################################################
# Sun MP Definitions
##############################################################################

ifeq ($(ARCH),SUNMP)
	X11LDFLAGS = -L/usr/local/X11/lib -R/usr/local/X11/lib
	ifeq ($(LP_SOLVER),CPLEX)
	   LPSOLVER_DEFS = -DSYSGNUSOLARIS
	endif
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM 
	SYSLIBS = -lsocket -lnsl
endif

##############################################################################
# X86 Solaris Definitions
##############################################################################

ifeq ($(ARCH),X86SOL2)
	X11LDFLAGS = -L/usr/local/X11/lib -R/usr/local/X11/lib
	ifeq ($(LP_SOLVER),CPLEX)
	   LPSOLVER_DEFS = -DSYSGNUSOLARIS
	endif
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM
	SYSLIBS = -lsocket -lnsl
endif

##############################################################################
# Alpha Definitions
##############################################################################

ifeq ($(ARCH),ALPHA)
	X11LDFLAGS = -L/usr/local/X11/lib
	MACH_DEP = -DHAS_RANDOM -DHAS_SRANDOM
	SYSLIBS =
endif

##############################################################################
##############################################################################
# !!!!!!!!!!!!!!!!!!!USER SHOULD NOT EDIT BELOW THIS LINE !!!!!!!!!!!!!!!!!!!!
##############################################################################
##############################################################################

##############################################################################
# Set the VER to "g" if using gcc
##############################################################################

ifeq ($(CC),gcc)
	VER = g
else
	VER = x
endif
ifeq ($(VER),g)
	VERSION=GNU
else
	VERSION=NOGNU
endif

##############################################################################
##############################################################################
# What to make ? This has to go here in case the user has any targets.
##############################################################################
##############################################################################

WHATTOMAKE = masterlib master
PWHATTOMAKE = pmaster
QWHATTOMAKE = pmaster
ifeq ($(COMPILE_IN_TM),FALSE)
WHATTOMAKE += tmlib tm
PWHATTOMAKE += ptm
QWHATTOMAKE += qtm
endif
ifeq ($(COMPILE_IN_LP),FALSE)
WHATTOMAKE += lplib lp
PWHATTOMAKE += plp
QWHATTOMAKE += qlp
endif
ifeq ($(COMPILE_IN_CP),FALSE)
WHATTOMAKE += cplib cp
PWHATTOMAKE += pcp
QWHATTOMAKE += qcp
endif
ifeq ($(COMPILE_IN_CG),FALSE)
WHATTOMAKE += cglib cg
PWHATTOMAKE += pcg
QWHATTOMAKE += qcg
endif

all : 
	$(MAKE) $(WHATTOMAKE)

pall :
	$(MAKE) $(PWHATOTOMAKE)

qall :
	$(MAKE) $(QWHATOTOMAKE)

##############################################################################
##############################################################################
#  Include the user specific makefile
##############################################################################
##############################################################################

include $(USERROOT)/Makefile

##############################################################################
##############################################################################
# Paths
##############################################################################
##############################################################################

##############################################################################
# Set the configuration path
##############################################################################

ifeq ($(COMPILE_IN_TM),TRUE)
	CONFIG:=1
else
	CONFIG:=0
endif 
ifeq ($(COMPILE_IN_LP),TRUE)
	CONFIG:=1$(CONFIG)
else
	CONFIG:=0$(CONFIG)
endif 
ifeq ($(COMPILE_IN_CG),TRUE)
	CONFIG:=1$(CONFIG)
else
	CONFIG:=0$(CONFIG)
endif 
ifeq ($(COMPILE_IN_CP),TRUE)
	CONFIG:=1$(CONFIG)
else
	CONFIG:=0$(CONFIG)
endif 

INCDIR  	 = $(EXTRAINCDIR) -I$(ROOT)/include -I$(USERROOT)/include 
INCDIR 		+= $(USER_INCDIR)
#__BEGIN_EXPERIMENTAL_SECTION__#
ifeq ($(DECOMP),TRUE)
INCDIR 		+= -I$(ROOT)/include/decomp
endif
#___END_EXPERIMENTAL_SECTION___#
OBJDIR		 = $(ROOT)/objects.$(ARCH)/$(CONFIG)
USER_OBJDIR  	 = $(USERROOT)/objects.$(ARCH)/$(CONFIG)
DEPDIR  	 = $(ROOT)/dep.$(ARCH)
USER_DEPDIR  	 = $(USERROOT)/dep.$(ARCH)
LIBDIR		 = $(ROOT)/lib.$(ARCH)
BINDIR  	 = $(USERROOT)/bin.$(ARCH)
SRCDIR  = \
	$(ROOT)/Common     : $(USERROOT)/Common    :\
	$(ROOT)/LP         : $(USERROOT)/LP        :\
	$(ROOT)/CutGen     : $(USERROOT)/CutGen    :\
	$(ROOT)/CutPool    : $(USERROOT)/CutPool   :\
	$(ROOT)/SolPool    : $(USERROOT)/SolPool   :\
	$(ROOT)/DrawGraph  : $(USERROOT)/DrawGraph :\
	$(ROOT)/Master     : $(USERROOT)/Master    :\
	$(ROOT)/include    : $(USERROOT)/include   :\
	$(ROOT)            : $(USERROOT)           :\
	$(ROOT)/TreeManager                        :\
	$(USER_SRC_PATH)

VPATH  = $(SRCDIR):$(USER_SRCDIR)

##############################################################################
# Put it together
##############################################################################

LDFLAGS = -L$(LIBDIR) $(X11LDFLAGS) $(COMMLDFLAGS) $(LPLDFLAGS) $(USERLDFLAGS)
EXTRAINCDIR = $(X11INCDIR) $(COMMINCDIR) $(LPINCDIR)

ifeq ($(CC),ompcc)
	LIBS  = -lX11 -lm -lompc -ltlog -lthread $(COMMLIBS) $(SYSLIBS) \
	$(USERLIBS)
else
	LIBS  = -lX11 -lm $(COMMLIBS) $(SYSLIBS) $(USERLIBS)
endif
ifeq ($(OPT),-O)
    ifeq ($(CC),gcc)
	OPT = -O3 
    endif
    ifeq ($(ARCH),RS6K)
	ifeq ($(CC),xlC)
	    OPT = -O3 -qmaxmem=16384 -qarch=pwr2 -qtune=pwr2s
            OPT += -bmaxdata:0x80000000 -bloadmap:main.map
	endif
    endif
endif

ifeq ($(OPT),-g)
   ifeq ($(VERSION),GNU)
      OPT = -g 
      EFENCE = -lefence
   endif
   ifeq ($(ARCH),RS6K)
      ifeq ($(CC),xlC)
         OPT = -bmaxdata:0x80000000 -bloadmap:main.map
         OPT += -bnso -bnodelcsect -bI:/lib/syscalls.exp
         EFENCE = -lefence
      endif
   endif
endif

##############################################################################
##############################################################################
# Purify related flags
##############################################################################
##############################################################################

ifeq ($(ARCH),SUN4SOL2)
	PURIFYCACHEDIR=$(HOME)/purify-quantify/cache/SUN4SOL2
	PUREBIN = /home/purify/purify-4.1-solaris2/purify
endif
ifeq ($(ARCH),X86SOL2)
	PURIFYCACHEDIR=$(HOME)/purify-quantify/cache/X86SOL2
	PUREBIN = /home/purify/purify-4.1-solaris2/purify
endif

PFLAGS = -cache-dir=$(PURIFYCACHEDIR) -chain-length=10 \
	 -user-path=$(USERROOT)/bin.$(ARCH) \
         #-log-file=$(USERROOT)/purelog_%v.%p \
         #-mail-to-user=$(USER) # -copy-fd-output-to-logfile=1,2
PURIFY = $(PUREBIN) $(PFLAGS)

##############################################################################
##############################################################################
# Quantify related flags
##############################################################################
##############################################################################

ifeq ($(ARCH),SUN4SOL2)
	QUANTIFYCACHEDIR=$(HOME)/purify-quantify/cache/SUN4SOL2
	QUANTIFYBIN = /opts/pure/quantify-2.1-solaris2/quantify
endif
ifeq ($(ARCH),X86SOL2)
	QUANTIFYCACHEDIR=$(HOME)/purify-quantify/cache/X86SOL2
	QUANTIFYBIN = /opts/pure/quantify-2.1-solaris2/quantify
endif
QFLAGS   = -cache-dir=$(QUANTIFYCACHEDIR) 
QFLAGS  += -user-path=$(ROOT)/$(USERROOT)/bin.$(ARCH)
QUANTIFY = $(QUANTIFYBIN) $(QFLAGS)

##############################################################################
##############################################################################
#  Extensions for filenames for various configurations
##############################################################################
##############################################################################

ifeq ($(COMPILE_IN_CG),TRUE)
LPEXT = _cg
endif
ifeq ($(COMPILE_IN_CP),TRUE)
CPEXT = _cp
endif
ifeq ($(COMPILE_IN_LP),TRUE)
TMEXT = _lp$(LPEXT)$(CPEXT)
TMLPLIB = $(LPLIB)
endif
ifeq ($(COMPILE_IN_TM),TRUE)
MASTEREXT = _tm$(TMEXT)
ifeq ($(COMPILE_IN_LP),TRUE)
MASTERLPLIB = $(LPLIB)
endif
endif

##############################################################################
##############################################################################
# Putting together DEF's, FLAGS
##############################################################################
##############################################################################

SYSDEFINES  = -D__$(LP_SOLVER)__ -D__$(COMM_PROTOCOL)__ $(LPSOLVER_DEFS) 
SYSDEFINES += $(MACH_DEP)

BB_DEFINES  = $(USER_BB_DEFINES)
ifeq ($(ROOT_NODE_ONLY),TRUE)
BB_DEFINES += -DROOT_NODE_ONLY
endif
ifeq ($(TRACE_PATH),TRUE)
BB_DEFINES += -DTRACE_PATH 
endif
ifeq ($(CHECK_CUT_VALIDITY),TRUE)
BB_DEFINES += -DCHECK_CUT_VALIDITY
endif
ifeq ($(CHECK_LP),TRUE)
BB_DEFINES += -DCOMPILE_CHECK_LP
endif

ifeq ($(COMPILE_IN_CG),TRUE)
BB_DEFINES += -DCOMPILE_IN_CG
endif
ifeq ($(COMPILE_IN_CP),TRUE)
BB_DEFINES += -DCOMPILE_IN_CP
endif
ifeq ($(COMPILE_IN_LP),TRUE)
BB_DEFINES+= -DCOMPILE_IN_LP
endif
ifeq ($(COMPILE_IN_TM), TRUE)
BB_DEFINES += -DCOMPILE_IN_TM
endif
ifeq ($(COMPILE_FRAC_BRANCHING),TRUE)
BB_DEFINES  += -DCOMPILE_FRAC_BRANCHING
endif
ifeq ($(DO_TESTS),TRUE)
BB_DEFINES  += -DDO_TESTS
endif
ifeq ($(BIG_PROBLEM),TRUE)
BB_DEFINES  += -DBIG_PROBLEM
endif
ifeq ($(TM_BASIS_TESTS),TRUE)
BB_DEFINES  += -DTM_BASIS_TESTS
endif
ifeq ($(EXACT_MACHINE_HANDLING),TRUE)
BB_DEFINES  += -DEXACT_MACHINE_HANDLING
endif
ifeq ($(STATISTICS),TRUE)
BB_DEFINES  += -DSTATISTICS
endif
ifeq ($(PSEUDO_COSTS),TRUE)
BB_DEFINES  += -DPSEUDO_COSTS
endif

#__BEGIN_EXPERIMENTAL_SECTION__#
##############################################################################
# DECOMP related stuff 
##############################################################################

ifeq ($(DECOMP),TRUE)
BB_DEFINES += -DCOMPILE_DECOMP
endif

#___END_EXPERIMENTAL_SECTION___#
##############################################################################
# Compiler flags
##############################################################################

STRICT_CHECKING = TRUE

DEFAULT_FLAGS = $(OPT) $(SYSDEFINES) $(BB_DEFINES) $(INCDIR)

MORECFLAGS =

ifeq ($(STRICT_CHECKING),TRUE)
ifeq ($(VERSION),GNU)
	MORECFLAGS = -ansi -pedantic -Wall -Wid-clash-81 -Wpointer-arith -Wwrite-strings -Wstrict-prototypes -Wmissing-prototypes -Wnested-externs -Winline -fnonnull-objects #-pipe
endif
else
MOREFLAGS = 
endif

CFLAGS = $(DEFAULT_FLAGS) $(MORECFLAGS) $(MOREFLAGS)

##############################################################################
##############################################################################
# Global source files
##############################################################################
##############################################################################

MASTER_SRC	= master.c master_wrapper.c master_io.c
DG_SRC		= draw_graph.c

ifeq ($(COMPILE_IN_TM), TRUE)
TM_SRC		= tm_func.c tm_proccomm.c
else
TM_SRC          = treemanager.c tm_func.c tm_proccomm.c
endif
ifeq ($(COMPILE_IN_LP),TRUE)
TM_SRC         += lp_solver.c lp_varfunc.c lp_rowfunc.c lp_genfunc.c
TM_SRC         += lp_proccomm.c lp_wrapper.c lp_free.c
ifeq ($(PSEUDO_COSTS),TRUE)
TM_SRC         += lp_pseudo_branch.c
else
TM_SRC         += lp_branch.c
endif
ifeq ($(COMPILE_IN_CG),TRUE)
TM_SRC         += cg_func.c 
#__BEGIN_EXPERIMENTAL_SECTION__#
ifeq ($(DECOMP),TRUE)
TM_SRC		+= decomp.c decomp_lp.c
endif
#___END_EXPERIMENTAL_SECTION___#
endif
endif
ifeq ($(COMPILE_IN_CP),TRUE)
TM_SRC	       += cp_proccomm.c cp_func.c
endif
ifeq ($(COMPILE_IN_TM),TRUE)
MASTER_SRC     += $(TM_SRC)
endif

LP_SRC		= lp_solver.c lp_varfunc.c lp_rowfunc.c lp_genfunc.c \
		  lp_proccomm.c lp_wrapper.c lp.c lp_free.c
ifeq ($(PSEUDO_COSTS),TRUE)
TM_SRC         += lp_pseudo_branch.c
else
TM_SRC         += lp_branch.c
endif
ifeq ($(COMPILE_IN_CG),TRUE)
LP_SRC         += cg_func.c 
#__BEGIN_EXPERIMENTAL_SECTION__#
ifeq ($(DECOMP),TRUE)
TM_SRC		+= decomp.c decomp_lp.c
endif
#___END_EXPERIMENTAL_SECTION___#
endif

CP_SRC		= cut_pool.c cp_proccomm.c cp_func.c

CG_SRC		= cut_gen.c cg_proccomm.c cg_func.c

QSORT_SRC	= qsortucb.c qsortucb_i.c qsortucb_ii.c qsortucb_id.c \
		  qsortucb_di.c qsortucb_ic.c
TIME_SRC	= timemeas.c
PROCCOMM_SRC	= proccomm.c
PACKCUT_SRC	= pack_cut.c
PACKARRAY_SRC	= pack_array.c

BB_SRC = $(MASTER_SRC) $(DG_SRC) $(TM_SRC) $(LP_SRC) $(CP_SRC) $(CG_SRC) \
	 $(QSORT_SRC) $(TIME_SRC) $(PROCCOMM_SRC) $(PACKCUT_SRC) \
	 $(PACKARRAY_SRC)

ALL_SRC = $(BB_SRC) $(USER_SRC)

##############################################################################
##############################################################################
# Global rules
##############################################################################
##############################################################################

$(OBJDIR)/%.o : %.c
	mkdir -p $(OBJDIR)
	@echo Compiling $*.c
	$(CC) $(CFLAGS) $(EFENCE_LD_OPTIONS) -c $< -o $@

$(DEPDIR)/%.d : %.c
	mkdir -p $(DEPDIR)
	@echo Creating dependency $*.d
	$(SHELL) -ec '$(CC) -MM $(CFLAGS) $< \
		| $(AWK) "(NR==1) {printf(\"$(OBJDIR)/$*.o \\\\\\n\"); \
                                 printf(\"$(DEPDIR)/$*.d :\\\\\\n \\\\\\n\"); \
                                } \
                        (NR!=1) {print;}" \
                > $@'

$(USER_OBJDIR)/%.o : %.c
	mkdir -p $(USER_OBJDIR)
	@echo Compiling $*.c
	$(CC) $(CFLAGS) $(EFENCE_LD_OPTIONS) -c $< -o $@

$(USER_DEPDIR)/%.d : %.c
	mkdir -p $(USER_DEPDIR)
	@echo Creating dependency $*.d
	$(SHELL) -ec '$(CC) -MM $(CFLAGS) $< \
		| $(AWK) "(NR==1) {printf(\"$(USER_OBJDIR)/$*.o \\\\\\n\"); \
                                 printf(\"$(USER_DEPDIR)/$*.d :\\\\\\n \\\\\\n\"); \
                                } \
                        (NR!=1) {print;}" \
                > $@'

##############################################################################
##############################################################################
# Master
##############################################################################
##############################################################################

ALL_MASTER	 = $(MASTER_SRC)
ALL_MASTER 	+= $(TIME_SRC)
ALL_MASTER 	+= $(QSORT_SRC)
ALL_MASTER 	+= $(PROCCOMM_SRC)
ALL_MASTER 	+= $(PACKCUT_SRC)
ALL_MASTER 	+= $(PACKARRAY_SRC)

MASTER_OBJS 	  = $(addprefix $(OBJDIR)/,$(notdir $(ALL_MASTER:.c=.o)))
MASTER_DEP 	  = $(addprefix $(DEPDIR)/,$(ALL_MASTER:.c=.d))
USER_MASTER_OBJS  = $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_MASTER_SRC:.c=.o)))
USER_MASTER_DEP   = $(addprefix $(USER_DEPDIR)/,$(USER_MASTER_SRC:.c=.d))

master : $(BINDIR)/master$(MASTEREXT)
	true

masterlib : $(LIBDIR)/libmaster$(MASTEREXT).a
	true

pmaster : $(BINDIR)/pmaster$(MASTEREXT)
	true

qmaster : $(BINDIR)/qmaster$(MASTEREXT)
	true

cmaster : $(BINDIR)/cmaster$(MASTEREXT)
	true

master_clean :
	cd $(OBJDIR)
	rm -f $(MASTER_OBJS)
	cd $(DEPDIR)
	rm -f $(MASTER_DEP)

master_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_MASTER_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_MASTER_DEP)

$(BINDIR)/master$(MASTEREXT) : $(USER_MASTER_DEP) $(USER_MASTER_OBJS) \
$(LIBDIR)/libmaster$(MASTEREXT).a 
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_MASTER_OBJS) \
	-lmaster$(MASTEREXT) $(LIBS) $(MASTERLPLIB) 
	@echo ""

$(LIBDIR)/libmaster$(MASTEREXT).a : $(MASTER_DEP) $(MASTER_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/libmaster$(MASTEREXT).a $(MASTER_OBJS)
	$(RANLIB) $(LIBDIR)/libmaster$(MASTEREXT).a
	@echo ""

$(BINDIR)/pmaster$(MASTEREXT) : $(USER_MASTER_DEP) $(USER_MASTER_OBJS) \
$(LIBDIR)/libmaster$(MASTEREXT).a 
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_MASTER_OBJS) \
	-lmaster$(MASTEREXT) $(MASTERLPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/qmaster$(MASTEREXT) : $(USER_MASTER_DEP) $(USER_MASTER_OBJS) \
$(LIBDIR)/libmaster$(MASTEREXT).a 
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(QUANTIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_MASTER_OBJS) \
	-lmaster$(MASTEREXT) $(MASTERLPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/cmaster$(MASTEREXT) : $(USER_MASTER_DEP) $(USER_MASTER_OBJS) \
$(LIBDIR)/libmaster$(MASTEREXT).a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CCMALLOC) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_MASTER_OBJS) \
	-lmaster$(MASTEREXT) $(MASTERLPLIB) $(LIBS)
	@echo ""

##############################################################################
##############################################################################
# DrawGraph
##############################################################################
##############################################################################

ALL_DG  = $(DG_SRC)
ALL_DG += $(PROCCOMM_SRC)

DG_OBJS 	= $(addprefix $(OBJDIR)/,$(notdir $(ALL_DG:.c=.o)))
DG_DEP  	= $(addprefix $(DEPDIR)/,$(ALL_DG:.c=.d))
USER_DG_OBJS 	= $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_DG_SRC:.c=.o)))
USER_DG_DEP  	= $(addprefix $(USER_DEPDIR)/,$(USER_DG_SRC:.c=.d))

dg : $(BINDIR)/dg
	true

dglib : $(LIBDIR)/libdg.a
	true

pdg : $(BINDIR)/pdg
	true

dg_clean :
	cd $(OBJDIR)
	rm -f $(DG_OBJS)
	cd $(DEPDIR)
	rm -f $(DG_DEP)

dg_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_DG_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_DG_DEP))

$(BINDIR)/dg : $(USER_DG_DEP) $(USER_DG_OBJS) $(LIBDIR)/libdg.a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_DG_OBJS) -ldg $(LIBS) 
	@echo ""

$(LIBDIR)/libdg.a : $(DG_DEP) $(DG_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/libdg.a $(DG_OBJS)
	$(RANLIB) $(LIBDIR)/libdg.a
	@echo ""

$(BINDIR)/pdg : $(USER_DG_DEP) $(USER_DG_OBJS) $(LIBDIR)/libdg.a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_DG_OBJS) -ldg \
	$(LIBS)
	@echo ""

##############################################################################
##############################################################################
# TreeManager
##############################################################################
##############################################################################

ALL_TM	 = $(TM_SRC)
ALL_TM 	+= $(TIME_SRC)
ALL_TM 	+= $(PROCCOMM_SRC)
ALL_TM 	+= $(PACKCUT_SRC)
ALL_TM 	+= $(PACKARRAY_SRC)
ifeq ($(COMPILE_IN_LP),TRUE)
ALL_TM  += $(QSORT_SRC)
endif

TM_OBJS 	= $(addprefix $(OBJDIR)/,$(notdir $(ALL_TM:.c=.o)))
TM_DEP  	= $(addprefix $(DEPDIR)/,$(ALL_TM:.c=.d))
USER_TM_OBJS 	= $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_TM_SRC:.c=.o)))
USER_TM_DEP  	= $(addprefix $(USER_DEPDIR)/,$(USER_TM_SRC:.c=.d))

tm : $(BINDIR)/tm$(TMEXT)
	true

tmlib : $(LIBDIR)/libtm$(TMEXT).a
	true

ptm : $(BINDIR)/ptm$(TMEXT)
	true

qtm : $(BINDIR)/qtm$(TMEXT)
	true

ctm : $(BINDIR)/ctm$(TMEXT)
	true

tm_clean :
	cd $(OBJDIR)
	rm -f $(TM_OBJS)
	cd $(DEPDIR)
	rm -f $(TM_DEP)

tm_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_TM_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_TM_DEP))

$(BINDIR)/tm$(TMEXT) : $(USER_TM_DEP) $(USER_TM_OBJS) $(LIBDIR)/libtm$(TMEXT).a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_TM_OBJS) -ltm$(TMEXT) \
	$(TMLPLIB) $(LIBS)
	@echo ""

$(LIBDIR)/libtm$(TMEXT).a : $(TM_DEP) $(TM_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/libtm$(TMEXT).a $(TM_OBJS)
	$(RANLIB) $(LIBDIR)/libtm$(TMEXT).a
	@echo ""

$(BINDIR)/ptm$(TMEXT) : $(USER_TM_DEP) $(USER_TM_OBJS) \
$(LIBDIR)/libtm$(TMEXT).a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_TM_OBJS) \
	-ltm$(TMEXT) $(TMLPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/qtm$(TMEXT) : $(USER_TM_DEP) $(USER_TM_OBJS) \
$(LIBDIR)/libtm$(TMEXT).a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(QUANTIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_TM_OBJS) \
	-ltm$(TMEXT) $(TMLPLIB) $(LIBS)
	@echo ""

$(BINDIR)/ctm$(TMEXT) : $(USER_TM_DEP) $(USER_TM_OBJS) \
$(LIBDIR)/libtm$(TMEXT).a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CCMALLOC) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_TM_OBJS) \
	-ltm$(TMEXT) $(TMLPLIB) $(LIBS)
	@echo ""

##############################################################################
##############################################################################
# LP
##############################################################################
##############################################################################

ALL_LP	 = $(LP_SRC)
ALL_LP 	+= $(TIME_SRC)
ALL_LP 	+= $(QSORT_SRC)
ALL_LP 	+= $(PROCCOMM_SRC)
ALL_LP 	+= $(PACKCUT_SRC)
ALL_LP 	+= $(PACKARRAY_SRC)

LP_OBJS 	= $(addprefix $(OBJDIR)/,$(notdir $(ALL_LP:.c=.o)))
LP_DEP 		= $(addprefix $(DEPDIR)/,$(ALL_LP:.c=.d))
USER_LP_OBJS 	= $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_LP_SRC:.c=.o)))
USER_LP_DEP 	= $(addprefix $(USER_DEPDIR)/,$(USER_LP_SRC:.c=.d))

lp : $(BINDIR)/lp$(LPEXT)
	true

lplib : $(LIBDIR)/liblp$(LPEXT).a
	true

plp : $(BINDIR)/plp$(LPEXT)
	true

qlp : $(BINDIR)/qlp$(LPEXT)
	true

clp : $(BINDIR)/clp$(LPEXT)
	true

lp_clean :
	cd $(OBJDIR)
	rm -f $(LP_OBJS)
	cd $(DEPDIR)
	rm -f $(LP_DEP)

lp_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_LP_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_LP_DEP))

$(BINDIR)/lp$(LPEXT) : $(USER_LP_DEP) $(USER_LP_OBJS) $(LIBDIR)/liblp$(LPEXT).a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_LP_OBJS) -llp$(LPEXT) \
	$(LPLIB) $(LIBS)
	@echo ""

$(LIBDIR)/liblp$(LPEXT).a : $(LP_DEP) $(LP_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/liblp$(LPEXT).a $(LP_OBJS)
	$(RANLIB) $(LIBDIR)/liblp$(LPEXT).a
	@echo ""

$(BINDIR)/plp$(LPEXT) : $(USER_LP_DEP) $(USER_LP_OBJS) \
$(LIBDIR)/liblp$(LPEXT).a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_LP_OBJS) \
	-llp$(LPEXT) $(LPLIB) $(LIBS)
	@echo ""

$(BINDIR)/qlp$(LPEXT) : $(USER_LP_DEP) $(USER_LP_OBJS) \
$(LIBDIR)/liblp$(LPEXT).a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(QUANTIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_LP_OBJS) \
	-llp$(LPEXT) $(LPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/clp$(LPEXT) : $(USER_LP_DEP) $(USER_LP_OBJS) \
$(LIBDIR)/liblp$(LPEXT).a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CCMALLOC) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_LP_OBJS) \
	-llp$(LPEXT) $(LPLIB) $(LIBS)
	@echo ""

##############################################################################
##############################################################################
# CutPool
##############################################################################
##############################################################################

ALL_CP	 = $(CP_SRC)
ALL_CP 	+= $(TIME_SRC)
ALL_CP 	+= $(QSORT_SRC)
ALL_CP 	+= $(PROCCOMM_SRC)
ALL_CP 	+= $(PACKCUT_SRC)

CP_OBJS 	= $(addprefix $(OBJDIR)/,$(notdir $(ALL_CP:.c=.o)))
CP_DEP  	= $(addprefix $(DEPDIR)/,$(ALL_CP:.c=.d))
USER_CP_OBJS 	= $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_CP_SRC:.c=.o)))
USER_CP_DEP  	= $(addprefix $(USER_DEPDIR)/,$(USER_CP_SRC:.c=.d))

cp : $(BINDIR)/cp
	true

cplib : $(LIBDIR)/libcp.a
	true

pcp : $(BINDIR)/pcp
	true

qcp : $(BINDIR)/qcp
	true

ccp : $(BINDIR)/ccp
	true

cp_clean :
	cd $(OBJDIR)
	rm -f $(CP_OBJS)
	cd $(DEPDIR)
	rm -f $(CP_DEP)

cp_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_CP_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_CP_DEP))

$(BINDIR)/cp : $(USER_CP_DEP) $(USER_CP_OBJS) $(LIBDIR)/libcp.a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CP_OBJS) -lcp $(LIBS) 
	@echo ""

$(LIBDIR)/libcp.a : $(CP_DEP) $(CP_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/libcp.a $(CP_OBJS)
	$(RANLIB) $(LIBDIR)/libcp.a
	@echo ""

$(BINDIR)/pcp : $(USER_CP_DEP) $(USER_CP_OBJS) $(LIBDIR)/libcp.a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CP_OBJS) -lcp \
	$(LIBS) 
	@echo ""

$(BINDIR)/qcp : $(USER_CP_DEP) $(USER_CP_OBJS) $(LIBDIR)/libcp.a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(QUANTIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CP_OBJS) -lcp \
	$(LIBS) 
	@echo ""

$(BINDIR)/ccp : $(USER_CP_DEP) $(USER_CP_OBJS) $(LIBDIR)/libcp.a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CCMALLOC) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CP_OBJS) -lcp \
	$(LIBS)
	@echo ""

##############################################################################
##############################################################################
# CutGen
##############################################################################
##############################################################################

ALL_CG	 = $(CG_SRC)
ALL_CG 	+= $(TIME_SRC)
ALL_CG 	+= $(QSORT_SRC)
ALL_CG 	+= $(PROCCOMM_SRC)
ALL_CG 	+= $(PACKCUT_SRC)

CG_OBJS = $(addprefix $(OBJDIR)/,$(notdir $(ALL_CG:.c=.o)))
CG_DEP  = $(addprefix $(DEPDIR)/,$(ALL_CG:.c=.d))
USER_CG_OBJS = $(addprefix $(USER_OBJDIR)/,$(notdir $(USER_CG_SRC:.c=.o)))
USER_CG_DEP  = $(addprefix $(USER_DEPDIR)/,$(USER_CG_SRC:.c=.d))

cg : $(BINDIR)/cg
	true

cglib : $(LIBDIR)/libcg.a
	true

pcg : $(BINDIR)/pcg
	true

qcg : $(BINDIR)/qcg
	true

ccg : $(BINDIR)/ccg
	true

cg_clean :
	cd $(OBJDIR)
	rm -f $(CG_OBJS)
	cd $(DEPDIR)
	rm -f $(CG_DEP)

cg_clean_user :
	cd $(USER_OBJDIR)
	rm -f $(USER_CG_OBJS)
	cd $(USER_DEPDIR)
	rm -f $(USER_CG_DEP))

$(BINDIR)/cg : $(USER_CG_DEP) $(USER_CG_OBJS) $(LIBDIR)/libcg.a
	@echo ""
	@echo "Linking $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CG_OBJS) -lcg $(LPLIB) $(LIBS) 
	@echo ""

$(LIBDIR)/libcg.a : $(CG_DEP) $(CG_OBJS)
	@echo ""
	@echo "Making $(notdir $@) ..."
	@echo ""
	mkdir -p $(LIBDIR)
	$(AR) $(LIBDIR)/libcg.a $(CG_OBJS)
	$(RANLIB) $(LIBDIR)/libcg.a
	@echo ""

$(BINDIR)/pcg : $(USER_CG_DEP) $(USER_CG_OBJS) $(LIBDIR)/libcg.a
	@echo ""
	@echo "Linking purified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(PURIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CG_OBJS) \
	-lcg $(LPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/qcg : $(USER_CG_DEP) $(USER_CG_OBJS) $(LIBDIR)/libcg.a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(QUANTIFY) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CG_OBJS) \
	-lcg $(LPLIB) $(LIBS) 
	@echo ""

$(BINDIR)/ccg : $(USER_CG_DEP) $(USER_CG_OBJS) $(LIBDIR)/libcg.a
	@echo ""
	@echo "Linking quantified $(notdir $@) ..."
	@echo ""
	mkdir -p $(BINDIR)
	$(CCMALLOC) $(CC) $(CFLAGS) $(LDFLAGS) -o $@ $(USER_CG_OBJS) \
	-lcg $(LPLIB) $(LIBS)
	@echo ""

###############################################################################
##############################################################################
# Special targets
##############################################################################
##############################################################################

.PHONY:	clean clean_all master_clean lp_clean cg_clean cp_clean tm_clean \
	dg_clean

clean :
	rm -rf $(ROOT)/objects.$(ARCH)

clean_user :
	rm -rf $(USERROOT)/objects.$(ARCH)

clean_dep :
	rm -rf $(DEPDIR)/ 

clean_user_dep :
	rm -rf $(USER_DEPDIR)/

clean_lib :
	rm -rf $(LIBDIR)

clean_bin :
	rm -rf $(BINDIR)

clean_all : clean clean_dep clean_user clean_user_dep clean_lib clean_bin
	true

.SILENT:
