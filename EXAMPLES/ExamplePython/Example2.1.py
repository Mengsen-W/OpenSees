# OpenSees -- Open System for Earthquake Engineering Simulation
# Pacific Earthquake Engineering Research Center
# http://opensees.berkeley.edu/
#
# Moment-Curvature Example 2.1
# ----------------------------
#  Zero length element with fiber section
#  Single Nodal Load, Static Analysis
# 
# Example Objectives
# ------------------
#  Moment-Curvature Analysis in OpenSees
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
import MomentCurvature

# ------------------------------
# Start of model generation
# ------------------------------

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
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

# Estimate yield curvature
# (Assuming no axial load and only top and bottom steel)
d = colDepth-cover;  # d -- from cover to rebar
epsy = fy/E;         # steel yield strain
Ky = epsy/(0.7*d)

# Print estimate to standard output
print("Estimated yield curvature: ", Ky)

# Set axial load 
P = -180.0

mu  = 15.0;     # Target ductility for analysis
numIncr = 100;  # Number of analysis increments

# Call the section analysis procedure
<<<<<<< HEAD
MomentCurvature.execute(ops, 1, P, Ky*mu, numIncr)

ops.wipe()
=======
MomentCurvature.execute(1, P, Ky*mu, numIncr)

wipe()
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
