# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Simply Supported Beam Example 6.1
# ---------------------------------
#  Simply supported beam modeled with
#  two dimensional solid elements
# 
# Example Objectives
# ------------------
#  test different quad elements
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

# create ModelBuilder (with two-dimensions and 2 DOF/node)
ops.model("BasicBuilder", "-ndm",2, "-ndf",2)

# set default units
ops.defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
=======
wipe()

# create ModelBuilder (with two-dimensions and 2 DOF/node)
model("BasicBuilder", "-ndm",2, "-ndf",2)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define the material
# -------------------
#                               matTag  E      nu      rho
<<<<<<< HEAD
ops.nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25, 6.75/g) 
=======
nDMaterial("ElasticIsotropic", 1, 1000.0, 0.25, 6.75/g) 
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define geometry
# ---------------
Quad = "quad"
#Quad = "SSPquad"
#Quad = "bbarQuad"
#Quad = "enhancedQuad"

thick = 2.0;
nx = 10; # NOTE: nx MUST BE EVEN FOR THIS EXAMPLE
ny = 4
bn = nx + 1 
l1 = int(nx/2 + 1)
l2 = int(l1 + ny*(nx+1))

# now create the nodes and elements using the block2D command
if (Quad == "quad" or Quad == "enhancedQuad"):
    #          numX numY startNode startEle eleType eleArgs? coords?
<<<<<<< HEAD
    ops.block2D(nx, ny, 1, 1,
=======
    block2D(nx, ny, 1, 1,
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
                Quad, thick, "PlaneStrain", 1,
                1,  0.0,  0.0,
                2, 40.0,  0.0,
                3, 40.0, 10.0,
                4,  0.0, 10.0)
elif (Quad == "SSPquad"):
    #          numX numY startNode startEle eleType eleArgs? coords?
<<<<<<< HEAD
    ops.block2D(nx, ny, 1, 1,
=======
    block2D(nx, ny, 1, 1,
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
                Quad, 1, "PlaneStrain", thick,
                1,  0.0,  0.0,
                2, 40.0,  0.0,
                3, 40.0, 10.0,
                4,  0.0, 10.0)
elif (Quad == "bbarQuad"):
    #          numX numY startNode startEle eleType eleArgs? coords?
<<<<<<< HEAD
    ops.block2D(nx, ny, 1, 1,
=======
    block2D(nx, ny, 1, 1,
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
                Quad, thick, 1,
                1,  0.0,  0.0,
                2, 40.0,  0.0,
                3, 40.0, 10.0,
                4,  0.0, 10.0)

# Single point constraints
#      node u1 u2    
<<<<<<< HEAD
ops.fix( 1, 1, 1)
ops.fix(bn, 0, 1)

# Define gravity loads
# create a Linear time series
ops.timeSeries("Linear", 1)
# create a Plain load pattern
ops.pattern("Plain", 1, 1, "-fact", 1.0)
ops.load(l1, 0.0, -1.0)
ops.load(l2, 0.0, -1.0)

# print model
#ops.printModel()
ops.printModel("-JSON", "-file", "Example6.1.json")
=======
fix( 1, 1, 1)
fix(bn, 0, 1)

# Define gravity loads
# create a Linear time series
timeSeries("Linear", 1)
# create a Plain load pattern
pattern("Plain", 1, 1, "-fact", 1.0)
load(l1, 0.0, -1.0)
load(l2, 0.0, -1.0)

# print model
#printModel()
printModel("-JSON", "-file", "Example6.1.json")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------- 
# End of model generation
# -----------------------


# --------------------------------------------------------------------
# Start of static analysis (creation of the analysis & analysis itself)
# --------------------------------------------------------------------

# create the system of equation
<<<<<<< HEAD
ops.system("ProfileSPD")

# create the DOF numberer
ops.numberer("RCM")

# create the constraint handler
ops.constraints("Plain")

# create the convergence test
ops.test("EnergyIncr", 1.0E-12, 10)

# create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm("Newton")

# create the load control with variable load steps
ops.integrator("LoadControl", 1.0, 1, 1.0, 10.0) 

# create the analysis object 
ops.analysis("Static")

# Perform the analysis
ops.analyze(10)     
=======
system("ProfileSPD")

# create the DOF numberer
numberer("RCM")

# create the constraint handler
constraints("Plain")

# create the convergence test
test("EnergyIncr", 1.0E-12, 10)

# create the solution algorithm, a Newton-Raphson algorithm
algorithm("Newton")

# create the load control with variable load steps
integrator("LoadControl", 1.0, 1, 1.0, 10.0) 

# create the analysis object 
analysis("Static")

# Perform the analysis
analyze(10)     
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# --------------------------
# End of static analysis
# --------------------------


# ----------------------------
# Start of recorder generation
# ----------------------------

<<<<<<< HEAD
ops.recorder("Node", "-file", "Node.out", "-time", "-node", l1, "-dof", 2, "disp")
#ops.recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#ops.recorder("display", g3, 10, 10, 800, 200, "-wipe")
=======
recorder("Node", "-file", "Node.out", "-time", "-node", l1, "-dof", 2, "disp")
#recorder("plot", "Node.out", "CenterNodeDisp", 625, 10, 625, 450, "-columns", 1, 2)

# create the display
#recorder("display", g3, 10, 10, 800, 200, "-wipe")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
#prp 20 5.0 100.0
#vup 0 1 0
#viewWindow -30 30 -10 10
#display 1 4 5

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

# create the system of equation
ops.system("ProfileSPD")

# create the DOF numberer
ops.numberer("RCM")

# create the constraint handler
ops.constraints("Plain")

# create the convergence test
ops.test("EnergyIncr", 1.0E-12, 10)

# create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm("Newton")

# create the integration scheme, the Newmark with gamma=0.5 and beta=0.25
ops.integrator("Newmark", 0.5, 0.25) 

# create the analysis object 
ops.analysis("Transient")

# record once at time 0
ops.record()

# Perform the transient analysis (20 sec)
#        numSteps  dt
ops.analyze(1000, 0.05)

ops.wipe()
=======
wipeAnalysis()
setTime(0.0)

# Now remove the loads and let the beam vibrate
remove("loadPattern", 1)

# create the system of equation
system("ProfileSPD")

# create the DOF numberer
numberer("RCM")

# create the constraint handler
constraints("Plain")

# create the convergence test
test("EnergyIncr", 1.0E-12, 10)

# create the solution algorithm, a Newton-Raphson algorithm
algorithm("Newton")

# create the integration scheme, the Newmark with gamma=0.5 and beta=0.25
integrator("Newmark", 0.5, 0.25) 

# create the analysis object 
analysis("Transient")

# record once at time 0
record()

# Perform the transient analysis (20 sec)
#        numSteps  dt
analyze(1000, 0.05)

wipe()
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
