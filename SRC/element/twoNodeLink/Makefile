
include ../../../Makefile.def
<<<<<<< HEAD
OBJS       = LinearElasticSpring.o \
    TwoNodeLink.o Inerter.o
=======
OBJS = Inerter.o \
       LinearElasticSpring.o \
	   TwoNodeLink.o
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Compilation control

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
