# Usage:
#   make            -- build with gfortran with all math optimisations
#   make FC=ifort   -- build with ifort with all math optimisations
#   make DEBUG=1    -- build with the debugging symbols, for gdb debugger
#   make FC=ifort DEBUG=1 -- build with debugging symbols using ifort.
#   To make and run regular model type make run, to make and run model where individuals 
#   are forced through an environment use make runex

# Supported Fortran compiler types
GF_FC = gfortran
GF_FC_W64 = x86_64-w64-mingw32-gfortran
IF_FC = ifort

# Choose the default compiler type, here it is gfortran. To use ifortan as
# the default set FC to $(IF_FC).
FC = $(GF_FC)

# Sources for the main program. Note that parameters.f90 go first in the
# list because the two remaining files / modules use the precision constant
# set in parameters.f90
SRC = parameters.f90 folder.f90 functions.f90 binary_IO.f90

# Main source file for the main program 
SRC_MAIN = Hormones_v2.f90

# Output executable
OUT = HORMONES.exe

# Source file for experimental program
SRC_EX = Hormones_v2_forcedE.f90

# Output executable
OUT_EX = HORMONES_EX.exe

# Options for GNU Fortran
# -Wall  does all checks in gfortran
GF_FFLAGS = -O3 -funroll-loops -fforce-addr -ffree-line-length-none -Wall -freal-4-real-8 -finteger-4-integer-8 -static-libgfortran -static -static-libgcc
# NOTE: the option that checks the array upper and lower bounds catches a
#       runtime error in the code:
#         At line 181 of file Hormones_v2.f90
#         Fortran runtime error: Index '0' of dimension 2 of array 'fitness'
#         below lower bound of 1
#       So it is disabled here:
#         -fbounds-check

# Options for Intel Fortran on Linux
#-fpe3 -warn -check bounds,pointers,format,uninit does checks in Intel Fortran
#-fpe3 on intel halts on arithmetic errors like NaN
IF_FFLAGS = -sox -O3 -ipo -fpe3 -warn -check bounds,pointers,format,uninit -autodouble

# On the Windows platform, intel fortran compiler uses a different pattern
# for compiler flags, starting with /Q
IF_FFLAGS_WINDOWS = /Qsox /O3 /Qipo /fpe3 /warn /check:bounds,pointers,format,uninit

# Check if debugging symbols should be included, then DEBUG turns off all
# optimisations and keeps debug symbols. The main option here is -g, it is
# the same in all compilers. To include debug symbols, make should be called
# with the DEBUG variable set to some value. For example:
#   make DEBUG=YES
#
# Note: -fbounds-check is absent/disabled as above to get around the BUG in
#       hormone model v2 code at line 181, fitness wrong array bounds.
ifdef DEBUG
	GF_FFLAGS = -O0 -g -ffree-line-length-none \
              -ffpe-trap=zero,invalid,overflow,underflow -Warray-temporaries \
              -Wall
	IF_FFLAGS = -O0 -g -fpe0 -warn -traceback -check bounds,pointers,format,uninit
	IF_FFLAGS_WINDOWS = /O0 /debug:all /fpe0 /warn /check:bounds,pointers,format,uninit
endif

#*******************************************************************************
# Determine what is the build platform, Windows / Unix
# ComSpec is defined on Windows, it may be case-sensitive, check with env.exe.
# Platform-specific commands for deleting and moving files are then defined.
ifdef ComSpec
	PLATFORM_TYPE=Windows
	RM ?= del
	MV ?= move
else
	PLATFORM_TYPE=Unix
	WHICH_CMD=which
	RM ?= rm -fr
	MV ?= mv -f
endif

# This is a special check for Intel fortran ifort running on Windows, on this
# platform combination, compiler flags have a different pattern.
ifeq ($(PLATFORM_TYPE)$(FC),Windowsifort)
	IF_FFLAGS:=$(IF_FFLAGS_WINDOWS)
endif

# Set build options depending on the specific compiler
ifeq ($(FC),$(GF_FC))
	FFLAGS = $(GF_FFLAGS)
endif

ifeq ($(FC),$(GF_FC_W64))
	FFLAGS = $(GF_FFLAGS)
endif

ifeq ($(FC),$(IF_FC))
	FFLAGS = $(IF_FFLAGS)
endif

# ==============================================================================
# Targets follow:
# Note: the 'all' and 'distclean' are the standard targets used in the
#       Code:Blocks project file supplied at the SVN. To get it into the
#       current directory with the code issue (man copy-paste):
#       @code
#        svn export https://svn.uib.no/aha-fortran/trunk/scripts/Project.cbp
#       @endcode
#       This project file uses this Makefile for building the program.

# This is a 'standard' default target.
all: $(OUT)

#This is the shortcut target take calls "make Hormones_v2_forcedE.f90"
ex: $(OUT_EX)

# This target produces the executable program.
$(OUT): $(SRC) $(SRC_MAIN)
	$(FC) $(FFLAGS) $^ -o $(OUT)

# This target produces the experimental executable program
$(OUT_EX): $(SRC) $(SRC_EX)
	$(FC) $(FFLAGS) $^ -o $(OUT_EX)

# Clean temporary mod object and output files, add more patterns as needed,
# e.g. output data files like *.bin
# Note: the dash before rm (-rm) says that the makefile should not stop if
#       the rm command finds any errors, such as cannot delete a file, just
#       keep going. It is okay for rm
.PHONY: clean
clean:
	-$(RM) $(OUT) $(OUT_EX) *.mod fort.* *.tmp *.o *.obj *.pdb *.ilk

# Distclean is also a 'standard' target. Here it just calls clean, and then
# additionally deletes binary output data files.
distclean: clean
	-$(RM) *bin

# Run the model, if the executable has 'outdated' status (the source has been
# edited since the last build), that will rebuild it first and then run.
run: $(OUT) $(SRC)
	$(info ****************************************************************)
	$(info *** Executing model $(OUT) now, runtime platform is $(PLATFORM_TYPE))
	$(info *** and compiler is $(FC))
	$(info ****************************************************************)
ifeq ($(PLATFORM_TYPE),Windows)
	$(OUT) $(RUNFLAG)
else
	./$(OUT) $(RUNFLAG)
endif
	$(info *** Model done)

# Run the experiment model, if the executable has 'outdated' status (the source has been
# edited since the last build), that will rebuild it first and then run.
runex: $(OUT_EX) $(SRC_EX)
	$(info ****************************************************************)
	$(info *** Executing model $(OUT_EX) now, runtime platform is $(PLATFORM_TYPE))
	$(info *** and compiler is $(FC))
	$(info ****************************************************************)
ifeq ($(PLATFORM_TYPE),Windows)
	$(OUT_EX) $(RUNFLAG)
else
	./$(OUT_EX) $(RUNFLAG)
endif
	$(info *** Model done)

ifeq ($(PLATFORM_TYPE),Windows)
	$(OUT_PVE) $(RUNFLAG)
else
	./$(OUT_PVE) $(RUNFLAG)
endif
	$(info *** Model done)
