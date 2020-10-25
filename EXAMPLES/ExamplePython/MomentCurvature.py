# A procedure for performing section analysis (only does
# moment-curvature, but can be easily modified to do any mode
# of section response.
#
# Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
# Date: June 2017
#
# Arguments
#	secTag -- tag identifying section to be analyzed
#	axialLoad -- axial load applied to section (negative is compression)
#	maxK -- maximum curvature reached during analysis
#	numIncr -- number of increments used to reach maxK (default 100)
#
# Sets up a recorder which writes moment-curvature results to file
# section$secTag.out ... the moment is in column 1, and curvature in column 2

<<<<<<< HEAD
def execute(ops, secTag, axialLoad, maxK, numIncr):
	
	# Define two nodes at (0,0)
	ops.node(1, 0.0, 0.0)
	ops.node(2, 0.0, 0.0)
	
	# Fix all degrees of freedom except axial and bending
	ops.fix(1, 1, 1, 1)
	ops.fix(2, 0, 1, 0)
	
	# Define element
	#                               tag ndI ndJ secTag
	ops.element("zeroLengthSection", 1, 1, 2, secTag)
	
	# Create recorder
	ops.recorder("Node", "-file", "section"+str(secTag)+".out", "-time", "-node", 2, "-dof", 3, "disp")
	
	# Define constant axial load
	ops.timeSeries("Constant", 1)
	ops.pattern("Plain", 1, 1)
	ops.load(2, axialLoad, 0.0, 0.0)
	
	# Define analysis parameters
	ops.system("BandGeneral")
	ops.numberer("Plain")
	ops.constraints("Plain")
	ops.test("NormUnbalance", 1.0E-9, 10)
	ops.algorithm("Newton")
	ops.integrator("LoadControl", 0.0)
	ops.analysis("Static")
	
	# Do one analysis for constant axial load
	ops.analyze(1)
	
	# Define reference moment
	ops.timeSeries("Linear", 2)
	ops.pattern("Plain", 2, 2)
	ops.load(2, 0.0, 0.0, 1.0)
=======
from opensees import *

def execute(secTag, axialLoad, maxK, numIncr):
	
	# Define two nodes at (0,0)
	node(1, 0.0, 0.0)
	node(2, 0.0, 0.0)
	
	# Fix all degrees of freedom except axial and bending
	fix(1, 1, 1, 1)
	fix(2, 0, 1, 0)
	
	# Define element
	#                               tag ndI ndJ secTag
	element("zeroLengthSection", 1, 1, 2, secTag)
	
	# Create recorder
	recorder("Node", "-file", "section"+str(secTag)+".out", "-time", "-node", 2, "-dof", 3, "disp")
	
	# Define constant axial load
	timeSeries("Constant", 1)
	pattern("Plain", 1, 1)
	load(2, axialLoad, 0.0, 0.0)
	
	# Define analysis parameters
	system("BandGeneral")
	numberer("Plain")
	constraints("Plain")
	test("NormUnbalance", 1.0E-9, 10)
	algorithm("Newton")
	integrator("LoadControl", 0.0)
	analysis("Static")
	
	# Do one analysis for constant axial load
	analyze(1)
	
	# Define reference moment
	timeSeries("Linear", 2)
	pattern("Plain", 2, 2)
	load(2, 0.0, 0.0, 1.0)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	
	# Compute curvature increment
	dK = maxK/numIncr
	
	# Use displacement control at node 2 for section analysis
<<<<<<< HEAD
	ops.integrator("DisplacementControl", 2, 3, dK, 1, dK, dK)
	
	# Do the section analysis
	ops.analyze(numIncr)
=======
	integrator("DisplacementControl", 2, 3, dK, 1, dK, dK)
	
	# Do the section analysis
	analyze(numIncr)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
	