# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Cantilever Beam Example 8.1
# ---------------------------
#  Cantilever beam modeled with
#  three dimensional brick elements
# 
# Example Objectives
# ------------------
#  test different brick elements
#  free vibration analysis starting from static deflection
#
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: September 2017

# import the OpenSees Python module
<<<<<<< HEAD
import opensees as ops
=======
from opensees import *
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------------
# Start of model generation
# ----------------------------

# remove existing model
<<<<<<< HEAD
ops.wipe()

# create ModelBuilder (with three-dimensions and 3 DOF/node)
ops.model("BasicBuilder", "-ndm",3, "-ndf",3)

# set default units
ops.defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
=======
wipe()

# create ModelBuilder (with three-dimensions and 3 DOF/node)
model("BasicBuilder", "-ndm",3, "-ndf",3)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define the material
# -------------------
#                               matTag  E     nu   rho
<<<<<<< HEAD
ops.nDMaterial("ElasticIsotropic", 1, 100.0, 0.25, 1.27) 
=======
nDMaterial("ElasticIsotropic", 1, 100.0, 0.25, 1.27) 
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define geometry
# ---------------
Brick = "stdBrick"
#Brick = "bbarBrick"
#Brick = "SSPbrick"

nz = 6
nx = 2 
ny = 2

nn = int((nz+1)*(nx+1)*(ny+1))

# mesh generation
#          numX numY numZ startNode startEle eleType eleArgs? coords?
<<<<<<< HEAD
ops.block3D(nx, ny, nz, 1, 1,
=======
block3D(nx, ny, nz, 1, 1,
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
            Brick, 1,
            1, -1.0, -1.0,  0.0,
            2,  1.0, -1.0,  0.0,
            3,  1.0,  1.0,  0.0,
            4, -1.0,  1.0,  0.0, 
            5, -1.0, -1.0, 10.0,
            6,  1.0, -1.0, 10.0,
            7,  1.0,  1.0, 10.0,
            8, -1.0,  1.0, 10.0)

# boundary conditions
<<<<<<< HEAD
ops.fixZ(0.0, 1, 1, 1) 

# Define point load
# create a Linear time series
ops.timeSeries("Linear", 1)
# create a Plain load pattern
load = 0.10
ops.pattern("Plain", 1, 1, "-fact", 1.0)
ops.load(nn, load, load, 0.0)

# print model
#ops.printModel()
ops.printModel("-JSON", "-file", "Example8.1.json")
=======
fixZ(0.0, 1, 1, 1) 

# Define point load
# create a Linear time series
timeSeries("Linear", 1)
# create a Plain load pattern
p = 0.10
pattern("Plain", 1, 1, "-fact", 1.0)
load(nn, p, p, 0.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example8.1.json")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------- 
# End of model generation
# -----------------------


# ------------------------
# Start of static analysis
# ------------------------

# Load control with variable load steps
#                            init  Jd  min   max
<<<<<<< HEAD
ops.integrator("LoadControl", 1.0, 1) 

# Convergence test
#                     tolerance maxIter displayCode
ops.test("NormUnbalance", 1.0E-10, 20, 0)

# Solution algorithm
ops.algorithm("Newton")

# DOF numberer
ops.numberer("RCM")

# Cosntraint handler
ops.constraints("Plain")

# System of equations solver
ops.system("ProfileSPD")

# Analysis for gravity load
ops.analysis("Static")

# Perform the analysis
ops.analyze(5)
=======
integrator("LoadControl", 1.0, 1) 

# Convergence test
#                     tolerance maxIter displayCode
test("NormUnbalance", 1.0E-10, 20, 0)

# Solution algorithm
algorithm("Newton")

# DOF numberer
numberer("RCM")

# Cosntraint handler
constraints("Plain")

# System of equations solver
system("ProfileSPD")

# Analysis for gravity load
analysis("Static")

# Perform the analysis
analyze(5)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# --------------------------
# End of static analysis
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

<<<<<<< HEAD
ops.recorder("Node", "-file", "Node.out", "-time", "-node", nn, "-dof", 1, "disp")
ops.recorder("Element", "-file", "Elem.out", "-time", "-eleRange", 1, 10, "material", "1", "strains")
#ops.recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#ops.recorder("display", "VibratingBeam", 100, 40, 500, 500, "-wipe")
=======
recorder("Node", "-file", "Node.out", "-time", "-node", nn, "-dof", 1, "disp")
recorder("Element", "-file", "Elem.out", "-time", "-eleRange", 1, 10, "material", "1", "strains")
#recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#recorder("display", "VibratingBeam", 100, 40, 500, 500, "-wipe")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
#prp -100 100 120.5
#vup 0 1 0 
#display 1 4 1 

# --------------------------
# End of recorder generation
# --------------------------


# ---------------------------------------
# Create and Perform the dynamic analysis
# ---------------------------------------

# Remove the static analysis & reset the time to 0.0
<<<<<<< HEAD
ops.wipeAnalysis()
ops.setTime(0.0)

# Now remove the loads and let the beam vibrate
ops.remove("loadPattern", 1)

# add some mass proportional damping
ops.rayleigh(0.01, 0.0, 0.0, 0.0)

# Create the transient analysis
ops.test("EnergyIncr", 1.0E-10, 20, 0)
ops.algorithm("Newton")
ops.numberer("RCM")
ops.constraints("Plain")
ops.system("ProfileSPD")
ops.integrator("Newmark", 0.5, 0.25)
ops.analysis("Transient")

# record once at time 0
ops.record()

# Perform the transient analysis (20 sec)
#         numSteps dt
ops.analyze(1000, 1.0)

ops.wipe()
=======
wipeAnalysis()
setTime(0.0)

# Now remove the loads and let the beam vibrate
remove("loadPattern", 1)

# add some mass proportional damping
rayleigh(0.01, 0.0, 0.0, 0.0)

# Create the transient analysis
test("EnergyIncr", 1.0E-10, 20, 0)
algorithm("Newton")
numberer("RCM")
constraints("Plain")
system("ProfileSPD")
integrator("Newmark", 0.5, 0.25)
analysis("Transient")

# record once at time 0
record()

# Perform the transient analysis (20 sec)
#         numSteps dt
analyze(1000, 1.0)

wipe()
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
