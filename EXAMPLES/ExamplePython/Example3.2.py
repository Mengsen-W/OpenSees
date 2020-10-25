# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Portal Frame Example 3.2
# ------------------------
#  Reinforced concrete one-bay, one-story frame
#  Distributed vertical load on girder
#  Lateral Load at top of frame
#  
# 
# Example Objectives
# -----------------
#  Nonlinear pushover analysis using Portal Frame Example 3.1 as starting point
# 
# Units: kips, in, sec
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: June 2017

# import the OpenSees Python module
<<<<<<< HEAD
import opensees as ops
=======
from opensees import *
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------------------------------------
# Start of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------

# remove existing model
<<<<<<< HEAD
ops.wipe()

# create ModelBuilder (with two-dimensions and 3 DOF/node)
ops.model("BasicBuilder", "-ndm",2, "-ndf",3)

# set default units
ops.defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
=======
wipe()

# create ModelBuilder (with two-dimensions and 3 DOF/node)
model("BasicBuilder", "-ndm",2, "-ndf",3)

# set default units
defaultUnits("-force", "kip", "-length", "in", "-time", "sec", "-temp", "F")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Create nodes
# ------------
# Set parameters for overall model geometry
width = 360.0
height = 144.0

# create nodes & add to Domain - command: node nodeId xCrd yCrd
<<<<<<< HEAD
ops.node(1, 0.0,   0.0)
ops.node(2, width, 0.0)
ops.node(3, 0.0,   height)
ops.node(4, width, height)

# set the boundary conditions - command: fix nodeID uxRestrnt? uyRestrnt? rzRestrnt?
ops.fix(1, 1, 1, 1)
ops.fix(2, 1, 1, 1)
=======
node(1, 0.0,   0.0)
node(2, width, 0.0)
node(3, 0.0,   height)
node(4, width, height)

# set the boundary conditions - command: fix nodeID uxRestrnt? uyRestrnt? rzRestrnt?
fix(1, 1, 1, 1)
fix(2, 1, 1, 1)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                        tag  f'c    ec0    f'cu   ecu
# Core concrete (confined)
<<<<<<< HEAD
ops.uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)
# Cover concrete (unconfined)
ops.uniaxialMaterial("Concrete01", 2, -5.0, -0.002, -0.0, -0.006)
=======
uniaxialMaterial("Concrete01", 1, -6.0, -0.004, -5.0, -0.014)
# Cover concrete (unconfined)
uniaxialMaterial("Concrete01", 2, -5.0, -0.002, -0.0, -0.006)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# STEEL
# Reinforcing steel 
fy = 60.0;      # Yield stress
E = 30000.0;    # Young's modulus
#                              tag fy  E0  b
<<<<<<< HEAD
ops.uniaxialMaterial("Steel01", 3, fy, E, 0.01)
=======
uniaxialMaterial("Steel01", 3, fy, E, 0.01)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define cross-section for nonlinear columns
# ------------------------------------------
# set some parameters
colWidth = 15.0
colDepth = 24.0 
cover = 1.5
As = 0.60;     # area of no. 7 bars

# some variables derived from the parameters
y1 = colDepth/2.0
z1 = colWidth/2.0

<<<<<<< HEAD
ops.section("Fiber", 1)
# Create the concrete core fibers
ops.patch("rect", 1, 10, 1, cover-y1, cover-z1, y1-cover, z1-cover)
# Create the concrete cover fibers (top, bottom, left, right)
ops.patch("rect", 2, 10, 1, -y1, z1-cover, y1, z1)
ops.patch("rect", 2, 10, 1, -y1, -z1, y1, cover-z1)
ops.patch("rect", 2,  2, 1, -y1, cover-z1, cover-y1, z1-cover)
ops.patch("rect", 2,  2, 1,  y1-cover, cover-z1, y1, z1-cover)
# Create the reinforcing fibers (left, middle, right)
ops.layer("straight", 3, 3, As, y1-cover, z1-cover, y1-cover, cover-z1)
ops.layer("straight", 3, 2, As, 0.0, z1-cover, 0.0, cover-z1)
ops.layer("straight", 3, 3, As, cover-y1, z1-cover, cover-y1, cover-z1)
# define beam integration
np = 5;  # number of integration points along length of element
ops.beamIntegration("Lobatto", 1, 1, np)
=======
section("Fiber", 1)
# Create the concrete core fibers
patch("rect", 1, 10, 1, cover-y1, cover-z1, y1-cover, z1-cover)
# Create the concrete cover fibers (top, bottom, left, right)
patch("rect", 2, 10, 1, -y1, z1-cover, y1, z1)
patch("rect", 2, 10, 1, -y1, -z1, y1, cover-z1)
patch("rect", 2,  2, 1, -y1, cover-z1, cover-y1, z1-cover)
patch("rect", 2,  2, 1,  y1-cover, cover-z1, y1, z1-cover)
# Create the reinforcing fibers (left, middle, right)
layer("straight", 3, 3, As, y1-cover, z1-cover, y1-cover, cover-z1)
layer("straight", 3, 2, As, 0.0, z1-cover, 0.0, cover-z1)
layer("straight", 3, 3, As, cover-y1, z1-cover, cover-y1, cover-z1)
# define beam integration
np = 5;  # number of integration points along length of element
beamIntegration("Lobatto", 1, 1, np)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define column elements
# ----------------------
# Geometry of column elements
#                       tag 
<<<<<<< HEAD
ops.geomTransf("PDelta", 1)
=======
geomTransf("PDelta", 1)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Create the coulumns using Beam-column elements
#                   tag ndI ndJ transfTag integrationTag
eleType = "forceBeamColumn"
<<<<<<< HEAD
ops.element(eleType, 1, 1, 3, 1, 1)
ops.element(eleType, 2, 2, 4, 1, 1)
=======
element(eleType, 1, 1, 3, 1, 1)
element(eleType, 2, 2, 4, 1, 1)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define beam element
# -----------------------------
# Geometry of column elements
#                tag 
<<<<<<< HEAD
ops.geomTransf("Linear", 2)

# Create the beam element
#                               tag ndI ndJ  A     E       Iz   transfTag
ops.element("elasticBeamColumn", 3, 3, 4, 360.0, 4030.0, 8640.0, 2)
=======
geomTransf("Linear", 2)

# Create the beam element
#                               tag ndI ndJ  A     E       Iz   transfTag
element("elasticBeamColumn", 3, 3, 4, 360.0, 4030.0, 8640.0, 2)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Define gravity loads
# --------------------
# Set a parameter for the axial load
P = 180.0;                # 10% of axial capacity of columns

# create a Linear TimeSeries (load factor varies linearly with time) - command: timeSeries Linear $tag
<<<<<<< HEAD
ops.timeSeries("Linear", 1)

# create a Plain load pattern - command: pattern Plain $tag $timeSeriesTag { $loads }
ops.pattern("Plain", 1, 1, "-fact", 1.0)

# create the nodal load - command: load nodeID xForce yForce zMoment
ops.load(3, 0.0, -P, 0.0)
ops.load(4, 0.0, -P, 0.0)
=======
timeSeries("Linear", 1)

# create a Plain load pattern - command: pattern Plain $tag $timeSeriesTag { $loads }
pattern("Plain", 1, 1, "-fact", 1.0)

# create the nodal load - command: load nodeID xForce yForce zMoment
load(3, 0.0, -P, 0.0)
load(4, 0.0, -P, 0.0)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ------------------------------
# End of model generation
# ------------------------------


# ------------------------------
# Start of analysis generation
# ------------------------------

# create the system of equation
<<<<<<< HEAD
ops.system("BandGeneral")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
ops.numberer("RCM")

# create the constraint handler, a Plain handler is used as homo constraints
ops.constraints("Plain")

# create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
ops.test("NormDispIncr", 1.0E-12, 10, 3)

# create the solution algorithm, a Newton-Raphson algorithm
ops.algorithm("Newton")

# create the integration scheme, the LoadControl scheme using steps of 0.1
ops.integrator("LoadControl", 0.1)

# create the analysis object 
ops.analysis("Static")
=======
system("BandGeneral")

# create the DOF numberer, the reverse Cuthill-McKee algorithm
numberer("RCM")

# create the constraint handler, a Plain handler is used as homo constraints
constraints("Plain")

# create the convergence test, the norm of the residual with a tolerance of 
# 1e-12 and a max number of iterations of 10
test("NormDispIncr", 1.0E-12, 10, 3)

# create the solution algorithm, a Newton-Raphson algorithm
algorithm("Newton")

# create the integration scheme, the LoadControl scheme using steps of 0.1
integrator("LoadControl", 0.1)

# create the analysis object 
analysis("Static")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ------------------------------
# End of analysis generation
# ------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# perform the gravity load analysis, requires 10 steps to reach the load level
<<<<<<< HEAD
ops.analyze(10)
=======
analyze(10)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

print("Gravity load analysis completed\n")

# Set the gravity loads to be constant & reset the time in the domain
<<<<<<< HEAD
ops.loadConst("-time", 0.0)
=======
loadConst("-time", 0.0)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------------------------------------
# End of Model Generation & Initial Gravity Analysis
# ----------------------------------------------------


# ----------------------------------------------------
# Start of additional modelling for lateral loads
# ----------------------------------------------------

# Define lateral loads
# --------------------
# Set some parameters
H = 10.0;		# Reference lateral load

# Set lateral load pattern with a Linear TimeSeries
<<<<<<< HEAD
ops.pattern("Plain", 2, 1, "-fact", 1.0)

# create the nodal load - command: load nodeID xForce yForce zMoment
ops.load(3, H, 0.0, 0.0)
ops.load(4, H, 0.0, 0.0)
=======
pattern("Plain", 2, 1, "-fact", 1.0)

# create the nodal load - command: load nodeID xForce yForce zMoment
load(3, H, 0.0, 0.0)
load(4, H, 0.0, 0.0)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------------------------------------
# End of additional modelling for lateral loads
# ----------------------------------------------------



# ----------------------------------------------------
# Start of modifications to analysis for push over
# ----------------------------------------------------

# Set some parameters
dU = 0.1;	        # Displacement increment

# Change the integration scheme to be displacement control
#                                    node dof init Jd min max
<<<<<<< HEAD
ops.integrator("DisplacementControl", 3, 1, dU, 1, dU, dU)
=======
integrator("DisplacementControl", 3, 1, dU, 1, dU, dU)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# ----------------------------------------------------
# End of modifications to analysis for push over
# ----------------------------------------------------


# ------------------------------
# Start of recorder generation
# ------------------------------

# Create a recorder to monitor nodal displacements
<<<<<<< HEAD
ops.recorder("Node", "-file", "node32.out", "-time", "-node", 3, 4, "-dof", 1, 2, 3, "disp")
#recorder plot node32.out hi 10 10 300 300 -columns 2 1

# Create a recorder to monitor element forces in columns
ops.recorder("EnvelopeElement", "-file", "ele32.out", "-time", "-ele", 1, 2, "localForce")
=======
recorder("Node", "-file", "node32.out", "-time", "-node", 3, 4, "-dof", 1, 2, 3, "disp")
#recorder plot node32.out hi 10 10 300 300 -columns 2 1

# Create a recorder to monitor element forces in columns
recorder("EnvelopeElement", "-file", "ele32.out", "-time", "-ele", 1, 2, "localForce")
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# --------------------------------
# End of recorder generation
# --------------------------------


# ------------------------------
# Finally perform the analysis
# ------------------------------

# record once at time 0
<<<<<<< HEAD
ops.record()
=======
record()
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Set some parameters
maxU = 15.0;	        # Max displacement
numSteps = int(maxU/dU)

# Perform the analysis
<<<<<<< HEAD
ok = ops.analyze(numSteps)

if (ok != 0):

    currentDisp = ops.nodeDisp(3, 1)
    ok = 0
    while ((ok == 0) and (currentDisp < maxU)):

        ok = ops.analyze(1)
=======
ok = analyze(numSteps)

if (ok != 0):

    currentDisp = nodeDisp(3, 1)
    ok = 0
    while ((ok == 0) and (currentDisp < maxU)):

        ok = analyze(1)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

        # if the analysis fails try initial tangent iteration
        if (ok != 0):
            print("regular newton failed .. lets try an initial stiffness for this step")
<<<<<<< HEAD
            ops.test("NormDispIncr", 1.0E-12, 1000) 
            ops.algorithm("ModifiedNewton", "-initial")
            ok = ops.analyze(1)
            if (ok == 0):
                print("that worked .. back to regular newton")
            ops.test("NormDispIncr", 1.0E-12, 10) 
            ops.algorithm("Newton")

        currentDisp = ops.nodeDisp(3, 1)
=======
            test("NormDispIncr", 1.0E-12, 1000) 
            algorithm("ModifiedNewton", "-initial")
            ok = analyze(1)
            if (ok == 0):
                print("that worked .. back to regular newton")
            test("NormDispIncr", 1.0E-12, 10) 
            algorithm("Newton")

        currentDisp = nodeDisp(3, 1)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
        
# Print a message to indicate if analysis successful or not
if (ok == 0):
    print("\nPushover analysis completed SUCCESSFULLY\n")
else:
    print("\nPushover analysis FAILED\n")

# Print the state at node 3
<<<<<<< HEAD
ops.printModel("node", 3)
ops.wipe()
=======
printModel("node", 3)
wipe()
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
