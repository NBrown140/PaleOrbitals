################################################################################
# Makefile for the Insolation Quantities
################################################################################
# fortran 90 compiler is required
# all fortran sources are in fixed form
########################################################################

#  Program macros
FFLAGS = 
FOR    = f90  -assume byterecl  # tru64 unix  - compaq fortran compiler
FOR    = ifort -fp-model precise # linux     - intel fortran compiler
#FOR    = xlf90  -qfixed        # AIX         - IBM fortran compiler
#FOR    = f90                   # macos X     - absoft fortran compiler
LIBS = -lm 

# Objects Files
OBJECTS= insolsub.o insola.o  prepinsol.o prepsub.o

.SUFFIXES: .f .o .par

# Compilation rules
.f.o: 
	$(FOR) $(FFLAGS) -c $<


############################################################
############################################################
# Actual Files
insola:	$(OBJECTS)
	$(FOR) $(FFLAGS)  $(OBJECTS) $(LIBS) -o insola

# Cleanup files 
clean:
	rm -f *.o *trace core *.mod work.pc work.pcl
	rm -f SOLCLI.BIN 

