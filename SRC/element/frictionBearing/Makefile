include ../../../Makefile.def

OBJS       = FlatSliderSimple2d.o \
	FlatSliderSimple3d.o \
	SingleFPSimple2d.o \
	SingleFPSimple3d.o \
	FPBearingPTV.o \
	TFP_Bearing.o \
	TFP_Bearing2d.o \
	TPB1D.o \
	TripleFrictionPendulum.o \
	MultiFP2d.o \
<<<<<<< HEAD
    RJWatsonEQS2d.o \
    RJWatsonEQS3d.o
=======
	RJWatsonEQS2d.o \
	RJWatsonEQS3d.o \
	TclFlatSliderCommand.o \
	TclRJWatsonEQSCommand.o \
	TclSingleFPCommand.o
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Compilation control

all:         $(OBJS)
	@$(CD) $(FE)/element/frictionBearing/frictionModel; $(MAKE);

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless
	@$(CD) $(FE)/element/frictionBearing/frictionModel; $(MAKE) wipe;

# DO NOT DELETE THIS LINE -- make depend depends on it.
