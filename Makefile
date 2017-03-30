#
# ---------------- Makefile for libUVLM: S.M. 06/09/2015 -----------------------
#-------------------------------------------------------------------------------
#
# Prerequisites:
# - update the path to eigen and boost libraries (EIGENDIR & BOOSTDIR variable).
#   Note that, in general, these do not require to be precompiled.
#
# Usage:
# - make all    : compile dynamic library and text executable
# - make so     : compile shared library only
# - make exe    : compile test executable only
# - make clean  : clean up
#
# Compiler:
# - compiler, compiling and linking options can be changed via the CPPCOMP,
#   LINKOPT and COMPOOPT. A summary of the option used with explainations is
#   provided
#
# Reference:
# http://oreilly.com/catalog/make3/book/
#
#-------------------------------------------------------------------------------

so: libUVLM

exe: mainexe

all: libUVLM mainexe

clean:
	rm -r -f ./lib/* ./bin/*


################################################################################
##  COMPILER OPTIONS.
#
#
# ---------------------------------------------------------------- Optimisation
# -O3: turn on all the optimisation options
# -Ofast: Disregard strict standards compliance. -Ofast enables all -O3 optimizations.
# It also enables optimizations that are not valid for all standard-compliant programs
#
# -------------------------------------------------------- Preprocessor Options:
# https://gcc.gnu.org/onlinedocs/gcc-3.0/gcc_3.html#SEC14
# 	-M: Instead of outputting the result of preprocessing, output a rule suitable
# for make describing the dependencies of the main source file. The preprocessor
# outputs one make rule containing the object file name for that source file, a
# colon, and the names of all the included files. Unless overridden explicitly,
# the object file name consists of the basename of the source file with any suffix
# replaced with object file suffix.
# 	-MD Like `-M' but the dependency information is written to a file rather than stdout.
# gcc will use the same file name and directory as the object file, but with the
# suffix `.d' instead.
#	-MMD: Like `-MD' except mention only user header files, not system -header files.
#	-MF: When used with `-M' or `-MM', specifies a file to write the dependencies to.
# This allows the preprocessor to write the preprocessed file to stdout normally.
# If no `-MF' switch is given, CPP sends the rules to stdout and suppresses normal
# preprocessed output.
# 	-MP: This option instructs CPP to add a phony target for each dependency
# other than the main file, causing each to depend on nothing. These dummy rules
# work around errors make gives if you remove header files without updating the
# Makefile to match.
# -MT By default CPP uses the main file name, including any path, and appends the
# object suffix, normally ".o", to it to obtain the name of the target for dependency
# generation. With `-MT' you can specify a target yourself, overriding the default one.
#
# -MMD, MP (kept)
# -MF, -MT (removed)
#
# -------------------------------------------------- Code Generation Conventions
# 	-fpic: Generate position-independent code (PIC) suitable for use in a shared library,
#if supported for the target machine. Such code accesses all constant addresses
#through a global offset table (GOT). The dynamic loader resolves the GOT entries
#when the program starts (the dynamic loader is not part of GCC; it is part of
#the operating system). If the GOT size for the linked executable exceeds a
#machine-specific maximum size, you get an error message from the linker
#indicating that `-fpic' does not work; in that case, recompile with `-fPIC'
#instead.
#	-fPIC: If supported for the target machine, emit position-independent code,
# suitable for dynamic linking and avoiding any limit on the size of the global
# offset table. This option makes a difference on the m68k, m88k, and the Sparc.
#
#
# ---------------------------------------------------------------------- Linking
# 	-shared: Produce a shared object which can then be linked with other objects
# to form an executable. Not all systems support this option. For predictable
# results, you must also specify the same set of options that were used to
# generate code (`-fpic', `-fPIC', or model suboptions) when you specify this option
# 	-fopenmp: activates the OpenMP extensions for C/C++ and Fortran. This enables
# the OpenMP directive #pragma omp in C/C++ and !$omp directives in free form,
# c$omp, *$omp and !$omp directives in fixed form, !$ conditional compilation
# sentinels in free form and c$, *$ and !$ sentinels in fixed form, for Fortran.
# The flag also arranges for automatic linking of the OpenMP runtime library
# (Runtime Library Routines).
#
#
# ------------------------------------------------------------------------ Other
# ref.: http://www.cs.fsu.edu/~jestes/howto/g++compiling.txt
#
# -c: translate source code into object code:
#		g++ -c xxx.cpp -----> xxx.o file
#
# -o: to name a target:
#		g++ -o yyy.o -c xxx.cpp ---> yyy.o file
#
# -Wall: enables a bunch of warning messages
#
# -g{level}: Request libging information and also use level to specify how much
# information.
# Level 1 produces minimal information, enough for making backtraces in parts of
# the program that you don't plan to lib. This includes descriptions of functions
# and external variables, but no information about local variables and no line numbers.
# Level 3 includes extra information, such as all the macro definitions present
# in the program. Some libgers support macro expansion when you use `-g3'.
#
#

CPPCOMP= g++-5

# include headers and external libraries
INCLOPT= -I"$(HDIR)" #-I"$(EIGENDIR)" -I"$(BOOSTDIR)"

FOMP=

# Linking options
LINKOPT= $(FOMP) -shared
# Compiling options
# try also -Ofast...
COMPOPT= -O0 -g3 -Wall -fmessage-length=0 $(FOMP) -fPIC -lm -std=gnu++14


################################################################################
## I/O SELECTION.
#
# Notes for Developer:
# - when inputting folder names, avoid spaces at the end of any assignment

# Building directories
BINDIR=./bin
LIBDIR=./lib
OBJDIR=$(LIBDIR)/obj

# Headers
HDIR=./include

# External Libraries
# EIGENDIR=/home/sm6110/git/eigen-eigen-3.2.5
# BOOSTDIR=/home/sm6110/git/boost_1_58_0

# Source Folder
SRCDIR=./src


################################################################################
## Source lists.
#
# This is required only to specify the prerequisites for the make targets

OBJS=\
	$(OBJDIR)/datatypesx.o \
	$(OBJDIR)/PanelTools.o \
	$(OBJDIR)/VLM.o \
	$(OBJDIR)/aicMats.o	\
	$(OBJDIR)/indices.o \
	$(OBJDIR)/triads.o \
	$(OBJDIR)/vorticity.o \
	$(OBJDIR)/wrapper.o


################################################################################
## Define make targets
#
# Automatic variables:
# $@: target filename
# $% filename element of an archive member specification
# $< filename of the first prerequisite
# $? names of all prerequisites that are newer than the target, separated by spaces.

#--------------------------------------------------------------- dynamic Library
all: mainexe

libUVLM: buildfolder $(OBJS)
	@echo 'Building target: $(LIBDIR)/libUVLM.so'
	$(CPPCOMP) $(LINKOPT) -o $(LIBDIR)/libUVLM.so $(OBJS)
	@echo 'Finished building target: $(BUILDDIR)/libUVLM.so'
	@echo 'Building target: ./bin/main.o'
	$(CPPCOMP) $(INCLOPT) $(COMPOPT) -c -o ./bin/main.o ./src/main.cpp
	@echo 'Finished building target: ./bin/main.o'

#	g++ -fopenmp -shared  -o ./lib/libUVLM.so ./lib/main.o  ./lib/lib/PanelTools.o ./lib/lib/VLM.o ./lib/lib/aicMats.o ./lib/lib/datatypesx.o ./lib/lib/indices.o ./lib/lib/triads.o ./lib/lib/vorticity.o ./lib/lib/wrapper.o
#	$(CPPCOMP)               $(LINKOPT) -o $(BUILDDIR)/libUVLM.so $(OBJS) ./lib/main.o

mainexe: buildfolder libUVLM $(OBJS)
	$(CPPCOMP) $(INCLOPT) $(COMPOPT) -c -o ./lib/main.o ./src/main.cpp
	$(CPPCOMP) $(INCLOPT) $(COMPOPT) -o $(BINDIR)/main ./lib/main.o $(OBJS)


################################################################################
## COMPILE DEPENDENCIES

# libraries
#$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	#$(CPPCOMP) $(INCLOPT) $(COMPOPT) -MF"$(OBJDIR)/$(%).d" -MT"$(OBJDIR)/$(%%).d" -o $@ $<
$(OBJDIR)/%.o: $(SRCDIR)/%.cpp
	@echo 'Building target: $@'
	$(CPPCOMP) $(INCLOPT) $(COMPOPT) -c $< -o $@
	@echo 'Finished building target: $@'
	@echo ' '

buildfolder:
	mkdir -p $(LIBDIR)
	mkdir -p $(BINDIR)
	mkdir -p $(OBJDIR)
#	@echo 'Source folder $(SRCDIR) contains:'
#	ls $(SRCDIR)
