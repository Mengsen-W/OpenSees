include ../../../../Makefile.def

OBJS       = BoundingCamClay.o \
             BoundingCamClay3D.o \
             BoundingCamClayPlaneStrain.o \
             ContactMaterial2D.o \
             ContactMaterial3D.o \
	     DruckerPrager.o \
             DruckerPrager3D.o \
	     DruckerPragerPlaneStrain.o \
	     InitialStateAnalysisWrapper.o\
	     ManzariDafalias.o\
	     ManzariDafalias3D.o\
	     ManzariDafaliasPlaneStrain.o\
	     ManzariDafaliasRO.o\
	     ManzariDafalias3DRO.o\
	     ManzariDafaliasPlaneStrainRO.o\
	     PM4Sand.o\
	     PM4Silt.o\
<<<<<<< HEAD
	     J2CyclicBoundingSurface.o
=======
	     J2CyclicBoundingSurface.o\
	     J2CyclicBoundingSurface3D.o\
	     J2CyclicBoundingSurfacePlaneStrain.o
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

all:	$(OBJS)

# Miscellaneous

tidy:
	@$(RM) $(RMFLAGS) Makefile.bak *~ #*# core

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.



