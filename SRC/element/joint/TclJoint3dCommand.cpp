/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
<<<<<<< HEAD
                                                                        
=======

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

// $Revision: 1.3 $
// $Date: 2004-09-01 04:01:27 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/TclJoint3dCommand.cpp,v $

// Written: Arash Altoontash, Gregory Deierlein,	Created: 04/01
// Revision: 
//				AAA		02/03
//
// Description: This file contains the implementation of the TclModelBuilder_addJoint3D()
// command. 
//
// What: "@(#) TclModelBuilder.C, revA"

#include <stdlib.h>
#include <string.h>
#include <Domain.h>

#include <Joint3D.h>
#include <TclModelBuilder.h>
#include <UniaxialMaterial.h>

//extern void printCommand(int argc, TCL_Char **argv);

int
<<<<<<< HEAD
TclModelBuilder_addJoint3D(ClientData clientData, Tcl_Interp *interp,  
			   int argc, 
			   TCL_Char **argv, 
			   Domain *theTclDomain,
			   TclModelBuilder *theTclBuilder)
=======
TclModelBuilder_addJoint3D(ClientData clientData, Tcl_Interp* interp,
  int argc,
  TCL_Char** argv,
  Domain* theTclDomain,
  TclModelBuilder* theTclBuilder)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
{
  // ensure the destructor has not been called
  if (theTclBuilder == 0) {
    opserr << "WARNING builder has been destroyed\n";
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
=======

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
  if (theTclBuilder->getNDM() != 3 || theTclBuilder->getNDF() != 6) {
    opserr << "WARNING -- model dimensions and/or nodal DOF not compatible with Joint3D element\n";
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  // check the number of arguments is correct
  int argStart = 2;
  
  if ( (argc-argStart) != 12 && (argc-argStart) != 16 ) {
=======

  // check the number of arguments is correct
  int argStart = 2;

  if ((argc - argStart) != 12 && (argc - argStart) != 16) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING incorrect number of arguments\n";
    //printCommand(argc, argv);
    opserr << "Want:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? MatX? MatY? MatZ? LrgDsp?\n";
    opserr << "or:\n";
    opserr << "element Joint3D Tag? NodI? NodJ? NodK? NodL? NodM? NodN? NodC? MatX? MatY? MatZ? LrgDsp? -damage DmgX DmgY DmgZ\n";
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
=======

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
  // get the id and end nodes
  int Joint3DId, iNode, jNode, kNode, lNode, mNode, nNode;
  if (Tcl_GetInt(interp, argv[argStart], &Joint3DId) != TCL_OK) {
    opserr << "WARNING invalid Joint3D eleTag" << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[1+argStart], &iNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[1 + argStart], &iNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid iNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[2+argStart], &jNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[2 + argStart], &jNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid jNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[3+argStart], &kNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[3 + argStart], &kNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid kNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[4+argStart], &lNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[4 + argStart], &lNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid lNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[5+argStart], &mNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[5 + argStart], &mNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid mNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  if (Tcl_GetInt(interp, argv[6+argStart], &nNode) != TCL_OK) {
=======

  if (Tcl_GetInt(interp, argv[6 + argStart], &nNode) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid nNode\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
  // Get the center node
  int CenterNodeTag;
<<<<<<< HEAD
  if (Tcl_GetInt(interp, argv[7+argStart], &CenterNodeTag) != TCL_OK) {
=======
  if (Tcl_GetInt(interp, argv[7 + argStart], &CenterNodeTag) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid tag for center node\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  // check domain for existence of internal node tag
  Node *CenterNode = theTclDomain->getNode(CenterNodeTag);
=======

  // check domain for existence of internal node tag
  Node* CenterNode = theTclDomain->getNode(CenterNodeTag);
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
  if (CenterNode != 0) {
    opserr << "WARNING node tag specified for the center node already exists.\n";
    opserr << "Use a new node tag.\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  UniaxialMaterial *MatX = NULL;
  int MatXid;
  if (Tcl_GetInt(interp, argv[8+argStart], &MatXid) != TCL_OK) {
=======

  UniaxialMaterial* MatX = NULL;
  int MatXid;
  if (Tcl_GetInt(interp, argv[8 + argStart], &MatXid) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid material ID for spring X\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  MatX = OPS_getUniaxialMaterial(MatXid);
  if ( MatX == NULL )
    {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << MatXid;
      opserr << "\nJoint3D element: " << Joint3DId << endln;
      return TCL_ERROR;
    }
  
  UniaxialMaterial *MatY = NULL;
  int MatYid;
  if (Tcl_GetInt(interp, argv[9+argStart], &MatYid) != TCL_OK) {
=======

  MatX = OPS_getUniaxialMaterial(MatXid);
  if (MatX == NULL)
  {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatXid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  UniaxialMaterial* MatY = NULL;
  int MatYid;
  if (Tcl_GetInt(interp, argv[9 + argStart], &MatYid) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid material ID for spring Y\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  MatY = OPS_getUniaxialMaterial(MatYid);
  if ( MatY == NULL )
    {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << MatYid;
      opserr << "\nJoint3D element: " << Joint3DId << endln;
      return TCL_ERROR;
    }		
  
  UniaxialMaterial *MatZ = NULL;
  int MatZid;
  if (Tcl_GetInt(interp, argv[10+argStart], &MatZid) != TCL_OK) {
=======

  MatY = OPS_getUniaxialMaterial(MatYid);
  if (MatY == NULL)
  {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatYid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  UniaxialMaterial* MatZ = NULL;
  int MatZid;
  if (Tcl_GetInt(interp, argv[10 + argStart], &MatZid) != TCL_OK) {
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    opserr << "WARNING invalid material ID for spring Z\n";
    opserr << "Joint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }
<<<<<<< HEAD
  
  MatZ = OPS_getUniaxialMaterial(MatZid);
  if ( MatZ == NULL )
    {
      opserr << "WARNING material not found\n";
      opserr << "Material: " << MatZid;
      opserr << "\nJoint3D element: " << Joint3DId << endln;
      return TCL_ERROR;
    }		
  
  int LargeDisp;
  if (Tcl_GetInt(interp, argv[11+argStart], &LargeDisp) != TCL_OK) {
    // use 0 as default
    LargeDisp = 0;
  }
  
  
  Joint3D *theJoint3D;
  // Decide to use which constructor, based on the number of arguments
  if ( (argc-argStart) == 12 ) {
    
    // Using Joint3D constructor without damage 
    theJoint3D = new Joint3D( Joint3DId,
			      iNode,jNode,kNode,lNode,mNode,nNode,CenterNodeTag,
			      *MatX,*MatY,*MatZ, theTclDomain, LargeDisp);
    
=======

  MatZ = OPS_getUniaxialMaterial(MatZid);
  if (MatZ == NULL)
  {
    opserr << "WARNING material not found\n";
    opserr << "Material: " << MatZid;
    opserr << "\nJoint3D element: " << Joint3DId << endln;
    return TCL_ERROR;
  }

  int LargeDisp;
  if (Tcl_GetInt(interp, argv[11 + argStart], &LargeDisp) != TCL_OK) {
    // use 0 as default
    LargeDisp = 0;
  }


  Joint3D* theJoint3D;
  // Decide to use which constructor, based on the number of arguments
  if ((argc - argStart) == 12) {

    // Using Joint3D constructor without damage 
    UniaxialMaterial* springModels[3] = { MatX, MatY, MatZ };
    theJoint3D = new Joint3D(Joint3DId,
      iNode, jNode, kNode, lNode, mNode, nNode, CenterNodeTag,
      springModels, theTclDomain, LargeDisp);

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    if (theJoint3D == 0) {
      opserr << "WARNING ran out of memory creating element\n";
      opserr << "Joint3D element: " << Joint3DId << endln;
      return TCL_ERROR;
    }
<<<<<<< HEAD
    
=======

>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    if (theTclDomain->addElement(theJoint3D) == false) {
      opserr << "WARNING could not add element to the domain\n";
      opserr << "Joint3D element: " << Joint3DId << endln;
      delete theJoint3D;
      return TCL_ERROR;
    }
<<<<<<< HEAD
    
    // if get here we have successfully created the element and added it to the domain
    return TCL_OK;
  }
  
  else 			// if ( (argc-argStart) == 16  )
    { 
      // Using Joint3D constructor with damage 
      // not implemented in this version
      return TCL_ERROR;
    }
=======

    // if get here we have successfully created the element and added it to the domain
    return TCL_OK;
  }

  else 			// if ( (argc-argStart) == 16  )
  {
    // Using Joint3D constructor with damage 
    // not implemented in this version
    return TCL_ERROR;
  }
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
  return TCL_ERROR;
}

