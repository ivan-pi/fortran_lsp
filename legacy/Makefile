#######################################################################
#                                                                     #
#  Usage:                                                             #
#                                                                     #
#  Type  make                  to compile all example programs        #
#  Type  make progN         to compile and run the example program #
#                              number N                               #
#  Type  make clean            to remove all but the generated output #
#                                                                     #
#######################################################################

#----------------------------------------------------------------------
# This section specifies the name of the linker and various
# compilation and linkage options.
#----------------------------------------------------------------------

FC=gfortran
FFLAGS= -Wall -fbacktrace -std=legacy
FLFLAGS=

PROGRAMS = prog1 prog2 prog3 prog4 prog5 prog6 prog7

#-------------------------------------------------------------------
# This section specifies which object files are needed for each
# example program.
#-------------------------------------------------------------------
prog1 = prog1.o gen.o h12.o
prog2 = prog2.o gen.o hfti.o h12.o diff.o
prog3 = prog3.o gen.o svdrs.o h12.o qrbd.o diff.o g1.o
prog4 = prog4.o sva.o svdrs.o h12.o qrbd.o diff.o g1.o mfeout.o
prog5 = prog5.o bndacc.o h12.o bndsol.o
prog6 = prog6.o svdrs.o h12.o qrbd.o diff.o g1.o ldp.o nnls.o
prog7 = prog7.o bvls.o

#--------------------------------------------------------------
# This is the default target.
#--------------------------------------------------------------

.PHONY : all
all : $(PROGRAMS)

#-------------------------------------------------------------
# This section shows how to link each example program.
#-------------------------------------------------------------
prog1 : $(prog1)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
prog2 : $(prog2)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
prog3 : $(prog3)
	$(FC) $(FFLAGS) -o $@ $(FLFLAGS) $^
prog4 : $(prog4)
	$(FC) $(FFLAGS)  -o $@ $(FLFLAGS) $^
prog5 : $(prog5)
	$(FC) $(FFLAGS)  -o $@ $(FLFLAGS) $^
prog6 : $(prog6)
	$(FC) $(FFLAGS)  -o $@ $(FLFLAGS) $^
prog7 : $(prog7)
	$(FC) $(FFLAGS)  -o $@ $(FLFLAGS) $^

#-------------------------------------------------------------
# This section contains default rule for creating Fortran
# object files.
#-------------------------------------------------------------

%.o: %.f
	$(FC) $(FFLAGS) -c $<
%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#----------------------------------------------
# This section shows how to clean up afterward.
#----------------------------------------------

.PHONY : clean cleanobj cleanmod
clean : cleanobj cleanmod
	rm -f $(PROGRAMS)
cleanobj :
	rm -f *.o
cleanmod :
	rm -f *.mod *__genmod.f90