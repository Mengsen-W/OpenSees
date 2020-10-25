include ../../../../Makefile.def

ifdef H5DRM_FLAG
<<<<<<< HEAD
  	H5_FILE = H5DRM.p
=======
	H5_FILE = H5DRM.o
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
else
	H5_FIle = 
endif

OBJS       = Mesh3DSubdomain.o \
	GeometricBrickDecorator.o \
	DRMBoundaryLayerDecorator.o \
	DRMLoadPattern.o \
	DRMLoadPatternWrapper.o \
	DRMInputHandler.o \
<<<<<<< HEAD
	PlaneDRMInputHandler.o $(H5_FILE)
=======
	PlaneDRMInputHandler.o \
	$(H5_FILE)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

all:         $(OBJS)

# Miscellaneous
tidy:	
	@$(RM) $(RMFLAGS) Makefile.bak *~ 

clean: tidy
	@$(RM) $(RMFLAGS) $(OBJS) *.o

spotless: clean

wipe: spotless

# DO NOT DELETE THIS LINE -- make depend depends on it.
