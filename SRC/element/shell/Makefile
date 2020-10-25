include ../../../Makefile.def

OBJS       = ShellMITC4.o \
	ShellMITC9.o \
	ShellDKGQ.o \
	ShellNLDKGQ.o \
	ShellDKGT.o \
	ShellNLDKGT.o \
	ShellMITC4Thermal.o \
	ShellNLDKGQThermal.o \
<<<<<<< HEAD
=======
	ASDShellQ4.o \
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	R3vectors.o \
	ShellANDeS.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
