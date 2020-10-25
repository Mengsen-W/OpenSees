include ../../../Makefile.def

OBJS       = FourNodeQuad.o \
	TclFourNodeQuadCommand.o \
	ConstantPressureVolumeQuad.o \
	EnhancedQuad.o \
	FourNodeQuad3d.o \
	NineNodeMixedQuad.o \
<<<<<<< HEAD
=======
	NineNodeQuad.o \
	EightNodeQuad.o \
	SixNodeTri.o \
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	FourNodeQuadWithSensitivity.o

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
