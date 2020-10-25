include ../../../../Makefile.def

OBJS       = Pinching.o \
	Bilinear.o \
<<<<<<< HEAD
	Bilin.o \
	Bilin02.o \
=======
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	Clough.o \
	CloughHenry.o \
	CloughDamage.o \
	PinchingDamage.o \
<<<<<<< HEAD
	Bilin.o \
=======
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	TclSnapMaterialCommand.o

all:         $(OBJS)

# Miscellaneous

tidy:   
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
