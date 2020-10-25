include ../../../Makefile.def

OBJS       = GenericClient.o \
<<<<<<< HEAD
	GenericCopy.o
=======
	GenericCopy.o \
	TclGenericClientCommand.o \
	TclGenericCopyCommand.o
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
