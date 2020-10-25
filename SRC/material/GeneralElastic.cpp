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

// Written: Mengsen Wang from Shenzhen University China
// Created: GuQuan by Xiamen University China
// Revision: user
//
// Description: This file contains the class definition for
// GeneralElastic.h adapted to user defined envelope
//   - Linear envelope with array
//
// What: "@(#) GeneralElastic.hh, revA"

#include <GeneralElastic.h>
#include <channel.h>
#include <elementAPI.h>
#include <tcl.h>

#include <cmath>

void *OPS_GeneralElasticMaterial(Tcl_Interp *interp, int argc,
                                 TCL_Char **argv) {
  if (argc < 3) {
    opserr << "WARNING insufficient arguments\n";
    opserr << "Want: uniaxialMaterial GeneralElasticMaterial tag? -strain {} "
              "-stress {}"
           << endln;
    return nullptr;
  }
  int tag;
  Vector *backboneStrain = nullptr;
  Vector *backboneStress = nullptr;

  if ((Tcl_GetInt(interp, argv[2], &tag) != 0)) {
    opserr << "Warring invaild tag in GeneralElastic Material" << endln;
    return nullptr;
  }

  if (strcmp(argv[3], "-strain") == 0) {
    int pathSize;
    TCL_Char **pathStrings;

    if (Tcl_SplitList(interp, argv[4], &pathSize, &pathStrings) != 0) {
      opserr << "Warning problem splitting path list " << argv[4] << " -";
      opserr << " stress -values {path} ...\n";
      return nullptr;
    }

    backboneStrain = new Vector(pathSize);

    for (int i = 0; i < pathSize; ++i) {
      double value;
      if (Tcl_GetDouble(interp, pathStrings[i], &value) != 0) {
        opserr << "Waning problem reading path data value " << pathStrings[i]
               << " - ";
        opserr << " -strain {} ...\n";
        return nullptr;
      }
      (*backboneStrain)(i) = value;
    }
  } else {
    opserr << "error command in generalElasticMaterial!" << endln;
    exit(-1);
  }

  if (strcmp(argv[5], "-stress") == 0) {
    int pathSize;
    TCL_Char **pathStrings;
    if (Tcl_SplitList(interp, argv[6], &pathSize, &pathStrings) != 0) {
      opserr << "Waring problem splitting path list " << argv[4] << " - ";
      opserr << " stress -values {path} ...\n";
      delete backboneStrain;
      return nullptr;
    }

    backboneStress = new Vector(pathSize);

    for (int i = 0; i < pathSize; ++i) {
      double value;
      if (Tcl_GetDouble(interp, pathStrings[i], &value) != 0) {
        opserr << "Waring problem reading path data value " << pathStrings[i]
               << " - ";
        opserr << " -strain {path} ...\n";

        delete backboneStrain;
        delete backboneStress;
        return nullptr;
      }
      (*backboneStress)(i) = value;
    }
  } else {
    opserr << "error command in generalElasticMaterial!" << endln;
    delete backboneStrain;
    exit(-1);
  }
  return static_cast<void *>(
      new GeneralElastic(tag, backboneStrain, backboneStress));
}

GeneralElastic::GeneralElastic(int tag, Vector *pBackboneStrain,
                               Vector *pbackboneStress)
    : UniaxialMaterial(tag, MAT_TAG_GeneralElastic),
      trialStrain(0.0),
      trialStress(0.0),
      tangent(0.0) {
  if (pBackboneStrain != nullptr) {
    backboneStrain = new Vector(*pBackboneStrain);
    backboneStress = new Vector(*pbackboneStress);
  } else {
    opserr << "Fatal: no backbone data in GeneralElastic Material!" << endln;
    exit(-1);
  }
  this->revertToStart();
}

GeneralElastic::GeneralElastic()
    : UniaxialMaterial(0, MAT_TAG_GeneralElastic),
      trialStrain(0.0),
      trialStress(0.0),
      tangent(0.0),
      trialStrainRate(0.0),
      backboneStrain(nullptr),
      backboneStress(nullptr) {
  this->revertToStart();
}

GeneralElastic::~GeneralElastic() {}

int GeneralElastic::setTrialStrain(double strain, double strainRate) {
  double sign = 1.0;
  if (fabs(strain) < 1.0e-14) sign = 0;
  if (strain < -1.0e-14) sign = -1;

  trialStrainRate = strain;
  trialStrainRate = strainRate;

  int i = 0;

  while ((i < backboneStrain->Size()) &&
         (fabs(trialStrain) > (*backboneStrain)[i] + 1.0e-14)) {
    ++i;
  }
  if (i == 0) {
    trialStress = 0.0;
    tangent = ((*backboneStress)[1] - (*backboneStress)[0]) /
              ((*backboneStrain)[1] - (*backboneStrain)[0]);
  } else if (i == backboneStrain->Size()) {
    // to big strain, keep horizontal
    trialStress = sign * (*backboneStress)[i - 1];
    tangent = 0.0;
  } else {
    // init stress + delta trial strain * delta(backboneStress) / delta(init
    // backboneStrain)

    tangent = ((*backboneStress)[i] - (*backboneStress)[i - 1]) /
              ((*backboneStrain)[i] - (*backboneStrain)[i - 1]);

    trialStress = (*backboneStress)[i - 1] +
                  (fabs(trialStrain) - (*backboneStrain)[i - 1]) * tangent;
    trialStress *= sign;
  }

  if (fabs(tangent) < 1.0e-14) {
    tangent = ((*backboneStress)[1] - (*backboneStress)[0]) /
              ((*backboneStrain)[1] - (*backboneStrain)[0]) * 1.0e-9;
  }
  if (tangent < -1.0e-14) {
    tangent = -((*backboneStress)[1] - (*backboneStress)[0]) /
              ((*backboneStrain)[1] - (*backboneStrain)[0]) * 1.0e-9;
  }

  return 0;
}

int GeneralElastic::setTrial(double strain, double &stress, double &theTangent,
                             double strainRate) {
  this->setTrialStrain(strain, strainRate);
  theTangent = tangent;
  stress = trialStress;
  return 0;
}

double GeneralElastic::getInitialTangent() {
  return ((*backboneStress)(1) - (*backboneStress)(0)) /
         ((*backboneStrain)(1) - (*backboneStrain)(0));
}

int GeneralElastic::revertToStart() {
  trialStrain = 0.0;
  trialStrainRate = 0.0;
  trialStress = 0.0;
  tangent = 0.0;

  return 0;
}

UniaxialMaterial *GeneralElastic::getCopy() {
  GeneralElastic *theCopy =
      new GeneralElastic(this->getTag(), backboneStrain, backboneStress);
  theCopy->trialStrain = trialStrain;
  theCopy->trialStrainRate = trialStrainRate;

  return theCopy;
}

void GeneralElastic::Print(OPS_Stream &s, int flag) {
  s << "Elastic tag: " << this->getTag() << endln;
  s << " backboneStrain: " << *backboneStrain << endln;
  s << " backboneStress: " << *backboneStress << endln;
}

Response *GeneralElastic::setResponse(const char **argv, int argc,
                                      OPS_Stream &theOutput) {
  Response *res = UniaxialMaterial::setResponse(argv, argc, theOutput);
  return res;
}

int GeneralElastic::getResponse(int responeseID, Information &matInfo) {
  return UniaxialMaterial::getResponse(responeseID, matInfo);
}

int GeneralElastic::sendSelf(int commitTag, Channel &theChannel) {
  int strain_size = backboneStrain->Size();
  int stress_size = backboneStress->Size();
  int basic = 7;

  static Vector data{basic + strain_size + stress_size};

  data(0) = this->getTag();
  data(1) = this->trialStrain;
  data(2) = this->trialStress;
  data(3) = this->tangent;
  data(4) = this->trialStrainRate;
  data(5) = strain_size;
  data(6) = stress_size;

  for (size_t i = 0; i < strain_size; ++i) {
    data(7 + i) = (*backboneStrain)(i);
  }

  for (size_t i = 0; i < stress_size; ++i) {
    data(7 + strain_size + i) = (*backboneStress)(i);
  }

  int res = res = theChannel.sendVector(this->getDbTag(), commitTag, data);

  if (res < 0) {
    opserr << "GerenalElsatic::sendSelf() - failed to send data\n";
  }

  return res;
}

int GeneralElastic::recvSelf(int commitTag, Channel &theChannel,
                             FEM_ObjectBroker &theBroker) {
  int strain_size = backboneStrain->Size();
  int stress_size = backboneStress->Size();
  int basic = 7;

  static Vector data{basic + strain_size + stress_size};

  int res = theChannel.recvVector(this->getDbTag(), commitTag, data);

  if (res < 0) {
    opserr << "GenernalElastic::recvSelf() - failed to receive data\n";
    return res;
  } else {
    this->setTag((int)data(0));
    this->trialStrain = data(1);
    this->trialStress = data(2);
    this->tangent = data(3);
    this->trialStrainRate = data(4);
    strain_size = data(5);
    stress_size = data(6);

    for (size_t i = 0; i < strain_size; ++i) {
      (*backboneStrain)(i) = data(7 + i);
    }

    for (size_t i = 0; i < stress_size; ++i) {
      (*backboneStress)(i) = data(7 + strain_size + i);
    }
  }

  return res;
}
