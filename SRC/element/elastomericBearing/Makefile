include ../../../Makefile.def

OBJS = ElastomericBearingPlasticity2d.o \
	ElastomericBearingPlasticity3d.o \
	ElastomericBearingBoucWen2d.o \
	ElastomericBearingBoucWen3d.o \
	ElastomericBearingBoucWenMod3d.o \
	ElastomericBearingUFRP2d.o \
	ElastomericX.o \
	HDR.o \
<<<<<<< HEAD
	LeadRubberX.o
=======
	LeadRubberX.o \
	TclElastomericBearingBoucWenCommand.o \
	TclElastomericBearingUFRPCommand.o \
	TclElastomericBearingPlasticityCommand.o
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
