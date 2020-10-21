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

// $Revision: 3.1.1 user version $
// $Date: 2020-10-20 21:49:52 $
// $Source: /OpenSees/SRC/material/uniaxial/GeneralElastic.h,v $

#ifndef GeneralElastic_h
#define GeneralElastic_h

#define MAT_TAG_GeneralElastic 20201020

#include <UniaxialMaterial.h>

// Written: Mengsen Wang from Shenzhen University China
// Created: GuQuan by Xiamen University China
// Revision: user
//
// Description: This file contains the class definition for
// GeneralElastic.h adapted to user defined envelope
//   - Linear envelope with array
//
// What: "@(#) GeneralElastic.hh, revA"

class GeneralElastic : public UniaxialMaterial {
 public:
  GeneralElastic(int tag, Vector *pBackboneStrain, Vector *pBackboneStress);
  GeneralElastic();
  ~GeneralElastic();

  const char *getClassType() const { return "GeneralElastic"; }

  int setTrialStrain(double strain, double strainRate = 0.0) override;
  int setTrial(double strain, double &stress, double &theTangent,
               double strainRate = 0.0);

  double getStrain() { return trialStrain; }
  double getStrainRate() { return trialStrainRate; }
  double getStress() { return trialStress; }
  double getTangent() { return tangent; }
  double getInitialTangent();
  int commitState() { return 0; }
  int revertToLastCommit() override { return 0; }
  int revertToStart();

  UniaxialMaterial *getCopy() override;

  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);

  void Print(OPS_Stream &s, int flag = 0);

  Response *setResponse(const char **argv, int argc, OPS_Stream &theOutput) override;
  int getResponse(int responseID, Information &matInfo) override;

 protected:
 private:
  double trialStrain;
  double trialStress;
  double tangent;
  double trialStrainRate;
  Vector *backboneStrain;
  Vector *backboneStress;
};

#endif  // GeneralElastic_h
