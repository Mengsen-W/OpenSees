include ../../../Makefile.def

OBJS       = DispBeamColumn2d.o \
	DispBeamColumn3d.o \
	DispBeamColumn2dThermal.o \
	DispBeamColumn3dThermal.o \
	DispBeamColumnNL2d.o \
	DispBeamColumn2dWithSensitivity.o \
	DispBeamColumn3dWithSensitivity.o \
	AxEqDispBeamColumn2d.o \
<<<<<<< HEAD
	TimoshenkoBeamColumn2d.o
=======
	TimoshenkoBeamColumn2d.o \
	DispBeamColumn3dID.o 
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0


all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
