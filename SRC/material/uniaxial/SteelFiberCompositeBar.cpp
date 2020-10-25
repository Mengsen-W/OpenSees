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
**   MengsenWang ShenZhenUniversity (mengsen_wang@163.com)            **
**                                                                    **
** ****************************************************************** */

// $Revision:  $3.1
// $Date: 2020-02-18 15:50:23 $
// $Source:
// /usr/local/cvs/OpenSees/SRC/material/uniaxial/SteelFiberCompositeBar.h,v $

// Written: MengsenWang
// Created: Feb 2020
//
// Description: This file contains the class definition for
// SteelFiberCompositeBar. SteelFiberCompositeBar provide the implementation
// of a one-dimensional Steel Fiber Composite Bar model with pinching of both
// force and deformation, damage due to deformation and energy, and degraded
// unkloading stiffness based on maximun ductility.

#include <Channel.h>
#include <OPS_Globals.h>
#include <SteelFiberCompositeBar.h>
#include <elementAPI.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>

void *OPS_SteelFiberCompositeBar(void) {
  // Pointer to a uniaxial material that will be return
  UniaxialMaterial *theMaterial = 0;

  int numArgs = OPS_GetNumRemainingInputArgs();
  if (numArgs != 12 && numArgs != 11) {
    opserr << "Want: uniaxialMaterial SteelFiberCompositeBar tag? mom1? rot1? "
              "mom2? rot2? mom3? rot3?"
           << "\npinchX? pinchY? damfc1? damfc2? <beta?> ";
    return 0;
  }

  // the position  of invalid number
  int iData[1];
  // all position number except tag
  double dData[11];
  for (int i = 0; i < 11; ++i) dData[i] = 0.0;

  int numData = 1;
  if (OPS_GetIntInput(&numData, iData) != 0) {
    opserr << "WARNING invalid tag for uniaxialMaterial SteelFiberCompositeBar"
           << endln;
    return 0;
  }

  //  remote tag
  numData = numArgs - 1;
  if (OPS_GetDoubleInput(&numData, dData) != 0) {
    opserr << "Invalid data for uniaxial SteelFiberCompositeBar" << iData[0]
           << endln;
    return 0;
  }

  // Parsing was successful, allocate the material
  theMaterial = new SteelFiberCompositeBar(
      iData[0], dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
      dData[6], dData[7], dData[8], dData[9], dData[10]);

  if (theMaterial == 0) {
    opserr << "WARNING could not creat uniaxialMaterial of type "
              "SteelFiberCompositeBar\n";
    return 0;
  }

  return theMaterial;
}

SteelFiberCompositeBar::SteelFiberCompositeBar(int tag, double mom1,
                                               double rot1, double mom2,
                                               double rot2, double mom3,
                                               double rot3, double pinchX,
                                               double pinchY, double damfc1,
                                               double damfc2, double beta)
    : UniaxialMaterial(tag, MAT_TAG_SteelFiberCompositeBar),
      pinchX(pinchX),
      pinchY(pinchY),
      damfc1(damfc1),
      damfc2(damfc2),
      beta(beta),
      mom1(mom1),
      rot1(rot1),
      mom2(mom2),
      rot2(rot2),
      mom3(mom3),
      rot3(rot3) {
  bool error = false;
  // Positive backbone parameters
  if (rot1 <= 0.0) error = true;
  if (rot2 <= 0.0) error = true;
  if (rot3 <= 0.0) error = true;
  if (error) {
    opserr << "SteelFiberCompositeBar::SteelFiberCompositeBar -- input "
              "backbone is not unique (one-to-one)\n";
    exit(-1);
  }

  energyA = 0.9 * (mom1 * rot1 + (mom1 + mom2) * (rot2 - rot1) +
                   (mom2 + mom3) * (rot3 - rot2));

  // Set envelope slopes
  this->setEnvelope();

  // Initialize history variables
  this->revertToStart();
  this->revertToLastCommit();
}

SteelFiberCompositeBar::SteelFiberCompositeBar()
    : UniaxialMaterial(0, MAT_TAG_SteelFiberCompositeBar),
      pinchX(0.0),
      pinchY(0.0),
      damfc1(0.0),
      damfc2(0.0),
      beta(0.0),
      mom1(0.0),
      rot1(0.0),
      mom2(0.0),
      rot2(0.0),
      mom3(0.0),
      rot3(0.0) {}

SteelFiberCompositeBar::~SteelFiberCompositeBar() {}

int SteelFiberCompositeBar::setTrialStrain(double strain, double strainRate) {
  if (TloadIndicator == 0 && strain == 0.0) return 0;
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TrotPu = CrotPu;
  TrotNu = CrotNu;

  Tstrain = strain;
  double dStrain = Tstrain - Cstrain;

  if (fabs(dStrain) < DBL_EPSILON) return 0;

  TloadIndicator = CloadIndicator;

  if (fabs(dStrain) < DBL_EPSILON) return 0;

  if (Tstrain >= CrotMax) {
    TrotMax = Tstrain;
    Ttangent = posEnvlpTangent(Tstrain);
    Tstress = posEnvlpStress(Tstrain);
    TloadIndicator = 1;
  } else {
    if (dStrain < 0.0)
      negativeIncrement(dStrain);
    else if (dStrain > 0.0)
      positiveIncrement(dStrain);
  }

  TenergyD = CenergyD + 0.5 * (Cstress + Tstress) * dStrain;

  return 0;
}

double SteelFiberCompositeBar::getStrain(void) { return Tstrain; }

double SteelFiberCompositeBar::getStress(void) { return Tstress; }

double SteelFiberCompositeBar::getTangent(void) { return Ttangent; }

void SteelFiberCompositeBar::positiveIncrement(double dStrain) {
  double kn = pow(CrotMin / rot1, beta);
  kn = (kn < 1.0) ? 1.0 : 1.0 / kn;
  double kp = pow(CrotMax / rot1, beta);
  kp = (kp < 1.0) ? 1.0 : 1.0 / kp;

  if (TloadIndicator == 2) {
    TloadIndicator = 1;
    if (Cstress <= 0.0) {
      TrotNu = Cstrain - Cstress / (Eu * kn);
      double energy = CenergyD - 0.5 * Cstress / (Eu * kn) * Cstress;
      double damfc = 0.0;
      if (CrotMin < rot1) {
        damfc = damfc2 * energy / energyA;
        damfc += damfc1 * (CrotMin - rot1) / rot1;
      }

      TrotMax = CrotMax * (1.0 + damfc);
    }
  }

  TloadIndicator = 1;

  if (TrotMax > POS_INF_STRAIN) TrotMax = POS_INF_STRAIN;

  TrotMax = (TrotMax > rot1) ? TrotMax : rot1;

  double maxmom = posEnvlpStress(TrotMax);
  double rotlim = negEnvlpRotlim(CrotMin);
  double rotrel = (rotlim > TrotNu) ? rotlim : TrotNu;

  double rotmp2 = TrotMax - (1.0 - pinchY) * maxmom / (Eu * kp);

  double rotch = rotrel + (rotmp2 - rotrel) * pinchX;

  double tmpmo1;
  double tmpmo2;

  if (Tstrain < TrotNu) {
    Ttangent = Eu * kn;
    Tstress = Cstress + Ttangent * dStrain;
    if (Tstress >= 0.0) {
      Tstress = 0.0;
      Ttangent = Eu * 1.0e-9;
    }
  }

  else if (Tstrain >= TrotNu && Tstrain < rotch) {
    if (Tstrain <= rotrel) {
      Tstress = 0.0;
      Ttangent = Eu * 1.0e-9;
    } else {
      Ttangent = maxmom * pinchY / (rotch - rotrel);
      tmpmo1 = Cstress + Eu * kp * dStrain;
      tmpmo2 = (Tstrain - rotrel) * Ttangent;
      if (tmpmo1 < tmpmo2) {
        Tstress = tmpmo1;
        Ttangent = Eu * kp;
      } else
        Tstress = tmpmo2;
    }
  }

  else {
    Ttangent = (1.0 - pinchY) * maxmom / (TrotMax - rotch);
    tmpmo1 = Cstress + Eu * kp * dStrain;
    tmpmo2 = pinchY * maxmom + (Tstrain - rotch) * Ttangent;
    if (tmpmo1 < tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = Eu * kp;
    } else
      Tstress = tmpmo2;
  }
}

void SteelFiberCompositeBar::negativeIncrement(double dStrain) {
  double kn = pow(CrotMin / rot1, beta);
  kn = (kn < 1.0) ? 1.0 : 1.0 / kn;
  double kp = pow(CrotMax / rot1, beta);
  kp = (kp < 1.0) ? 1.0 : 1.0 / kp;

  if (TloadIndicator == 1) {
    TloadIndicator = 2;
    if (Cstress >= 0.0) {
      TrotPu = Cstrain - Cstress / (Eu * kp);
      double energy = CenergyD - 0.5 * Cstress / (Eu * kp) * Cstress;
      double damfc = 0.0;
      if (CrotMax > rot1) {
        damfc = damfc2 * energy / energyA;
        damfc += damfc1 * (CrotMax - rot1) / rot1;
      }

      TrotMin = CrotMin * (1.0 + damfc);
    }
  }

  TloadIndicator = 2;

  if (TrotMin < NEG_INF_STRAIN) TrotMin = NEG_INF_STRAIN;

  TrotMin = (TrotMin < rot1) ? TrotMin : rot1;

  double minmom = negEnvlpStress(TrotMin);
  double rotlim = posEnvlpRotlim(CrotMax);
  double rotrel = (rotlim < TrotPu) ? rotlim : TrotPu;

  double rotmp2 = TrotMin - (1.0 - pinchY) * minmom / (Eu * kn);

  double rotch = rotrel + (rotmp2 - rotrel) * pinchX;

  double tmpmo1;
  double tmpmo2;

  if (Tstrain > TrotPu) {
    Ttangent = Eu * kp;
    Tstress = Cstress + Ttangent * dStrain;
    if (Tstress <= 0.0) {
      Tstress = 0.0;
      Ttangent = Eu * 1.0e-9;
    }
  }

  else if (Tstrain <= TrotPu && Tstrain > rotch) {
    if (Tstrain >= rotrel) {
      Tstress = 0.0;
      Ttangent = Eu * 1.0e-9;
    } else {
      Ttangent = minmom * pinchY / (rotch - rotrel);
      tmpmo1 = Cstress + Eu * kn * dStrain;
      tmpmo2 = (Tstrain - rotrel) * Ttangent;
      if (tmpmo1 > tmpmo2) {
        Tstress = tmpmo1;
        Ttangent = Eu * kn;
      } else
        Tstress = tmpmo2;
    }
  }

  else {
    Ttangent = (1.0 - pinchY) * minmom / (TrotMin - rotch);
    tmpmo1 = Cstress + Eu * kn * dStrain;
    tmpmo2 = pinchY * minmom + (Tstrain - rotch) * Ttangent;
    if (tmpmo1 > tmpmo2) {
      Tstress = tmpmo1;
      Ttangent = Eu * kn;
    } else
      Tstress = tmpmo2;
  }
}

int SteelFiberCompositeBar::commitState(void) {
  CrotMax = TrotMax;
  CrotMin = TrotMin;
  CrotPu = TrotPu;
  CrotNu = TrotNu;
  CenergyD = TenergyD;
  CloadIndicator = TloadIndicator;

  Cstress = Tstress;
  Cstrain = Tstrain;
  return 0;
}

int SteelFiberCompositeBar::revertToLastCommit(void) {
  TrotMax = CrotMax;
  TrotMin = CrotMin;
  TrotPu = CrotPu;
  TrotNu = CrotNu;
  TenergyD = CenergyD;
  TloadIndicator = CloadIndicator;

  Tstress = Cstress;
  Tstrain = Cstrain;

  return 0;
}

int SteelFiberCompositeBar::revertToStart(void) {
  CrotMax = 0.0;
  CrotMin = 0.0;
  CrotPu = 0.0;
  CrotNu = 0.0;
  CenergyD = 0.0;
  CloadIndicator = 0;

  Cstress = 0.0;
  Cstrain = 0.0;

  Tstrain = 0;
  Tstress = 0;
  Ttangent = E1;

  return 0;
}

UniaxialMaterial *SteelFiberCompositeBar::getCopy(void) {
  SteelFiberCompositeBar *theCopy = new SteelFiberCompositeBar(
      this->getTag(), mom1, rot1, mom2, rot2, mom3, rot3,
       pinchX, pinchY, damfc1, damfc2, beta);

  theCopy->CrotMax = CrotMax;
  theCopy->CrotMin = CrotMin;
  theCopy->CrotPu = CrotPu;
  theCopy->CrotNu = CrotNu;
  theCopy->CenergyD = CenergyD;
  theCopy->CloadIndicator = CloadIndicator;
  theCopy->Cstress = Cstress;
  theCopy->Cstrain = Cstrain;
  theCopy->Ttangent = Ttangent;

  return theCopy;
}

int SteelFiberCompositeBar::sendSelf(int commitTag, Channel &theChannel) {
  int res = 0;

  static Vector data(21);

  data(0) = this->getTag();
  data(1) = mom1;
  data(2) = rot1;
  data(3) = mom2;
  data(4) = rot2;
  data(5) = mom3;
  data(6) = rot3;
  data(7) = pinchX;
  data(8) = pinchY;
  data(9) = damfc1;
  data(10) = damfc2;
  data(11) = beta;
  data(12) = CrotMax;
  data(13) = CrotMin;
  data(14) = CrotPu;
  data(15) = CrotNu;
  data(16) = CenergyD;
  data(17) = CloadIndicator;
  data(18) = Cstress;
  data(19) = Cstrain;
  data(20) = Ttangent;

  res = theChannel.sendVector(this->getDbTag(), commitTag, data);
  if (res < 0)
    opserr << "SteelFiberCompositeBar::sendSelf() - failed to send data\n";

  return res;
}

int SteelFiberCompositeBar::recvSelf(int commitTag, Channel &theChannel,
                                 FEM_ObjectBroker &theBroker) {
  int res = 0;

  static Vector data(21);
  res = theChannel.recvVector(this->getDbTag(), commitTag, data);

  if (res < 0) {
    opserr << "SteelFiberCompositeBar::recvSelf() - failed to receive data\n";
    return res;
  } else {
    this->setTag((int)data(0));
    mom1 = data(1);
    rot1 = data(2);
    mom2 = data(3);
    rot2 = data(4);
    mom3 = data(5);
    rot3 = data(6);
    pinchX = data(7);
    pinchY = data(8);
    damfc1 = data(9);
    damfc2 = data(10);
    beta = data(11);

    CrotMax = data(12);
    CrotMin = data(13);
    CrotPu = data(14);
    CrotNu = data(15);
    CenergyD = data(16);
    CloadIndicator = int(data(17));
    Cstress = data(18);
    Cstrain = data(19);
    Ttangent = data(20);

    // set the trial values
    TrotMax = CrotMax;
    TrotMin = CrotMin;
    TrotPu = CrotPu;
    TrotNu = CrotNu;
    TenergyD = CenergyD;
    TloadIndicator = CloadIndicator;
    Tstress = Cstress;
    Tstrain = Cstrain;
  }

  // Set envelope slopes
  this->setEnvelope();

  return 0;
}

void SteelFiberCompositeBar::Print(OPS_Stream &s, int flag) {
  if (flag == OPS_PRINT_PRINTMODEL_MATERIAL) {
    s << "HSteelFiberCompositeBar, tag: " << this->getTag() << endln;
    s << "s1: " << mom1 << endln;
    s << "e1: " << rot1 << endln;
    s << "E1: " << E1 << endln;
    s << "s2: " << mom2 << endln;
    s << "e2: " << rot2 << endln;
    s << "E2: " << E2 << endln;
    s << "s3: " << mom3 << endln;
    s << "e3: " << rot3 << endln;
    s << "E3: " << E3 << endln;

    s << "pinchX: " << pinchX << endln;
    s << "pinchY: " << pinchY << endln;
    s << "damfc1: " << damfc1 << endln;
    s << "damfc2: " << damfc2 << endln;
    s << "energyA: " << energyA << endln;
    s << "beta: " << beta << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": \"" << this->getTag() << "\", ";
    s << "\"type\": \"SteelFiberCompositeBar\", ";
    s << "\"s1\": " << mom1 << ", ";
    s << "\"e1\": " << rot1 << ", ";
    s << "\"E1\": " << E1 << ", ";
    s << "\"s2\": " << mom2 << ", ";
    s << "\"e2\": " << rot2 << ", ";
    s << "\"E2\": " << E2 << ", ";
    s << "\"s3\": " << mom3 << ", ";
    s << "\"e3\": " << rot3 << ", ";
    s << "\"E3\": " << E3 << ", ";

    s << "\"pinchX\": " << pinchX << ", ";
    s << "\"pinchY\": " << pinchY << ", ";
    s << "\"damfc1\": " << damfc1 << ", ";
    s << "\"damfc2\": " << damfc2 << ", ";
    s << "\"energyA\": " << energyA << ", ";
    s << "\"beta\": " << beta << "}";
  }
}

void SteelFiberCompositeBar::setEnvelope(void) {
  E1 = mom1 / rot1;
  E2 = (mom2 - mom1) / (rot2 - rot1);
  E3 = (mom3 - mom2) / (rot3 - rot2);

  Eu = E1;
  if (E2 > Eu) Eu = E2;
  if (E3 > Eu) Eu = E3;

  Eu = E1;
  if (E2 > Eu) Eu = E2;
  if (E3 > Eu) Eu = E3;
}

double SteelFiberCompositeBar::posEnvlpStress(double strain) {
  if (strain <= 0.0)
    return 0.0;
  else if (strain <= rot1)
    return E1 * strain;
  else if (strain <= rot2)
    return mom1 + E2 * (strain - rot1);
  else if (strain <= rot3 || E3 > 0.0)
    return mom2 + E3 * (strain - rot2);
  else
    return mom3;
}

double SteelFiberCompositeBar::negEnvlpStress(double strain) {
  if (strain >= 0.0)
    return 0.0;
  else if (strain >= rot1)
    return E1 * strain;
  else if (strain >= rot2)
    return mom1 + E2 * (strain - rot1);
  else if (strain >= rot3 || E3 > 0.0)
    return mom2 + E3 * (strain - rot2);
  else
    return mom3;
}

double SteelFiberCompositeBar::posEnvlpTangent(double strain) {
  if (strain < 0.0)
    return E1 * 1.0e-9;
  else if (strain <= rot1)
    return E1;
  else if (strain <= rot2)
    return E2;
  else if (strain <= rot3 || E3 > 0.0)
    return E3;
  else
    return E1 * 1.0e-9;
}

double SteelFiberCompositeBar::negEnvlpTangent(double strain) {
  if (strain > 0.0)
    return E1 * 1.0e-9;
  else if (strain >= rot1)
    return E1;
  else if (strain >= rot2)
    return E2;
  else if (strain >= rot3 || E3 > 0.0)
    return E3;
  else
    return E1 * 1.0e-9;
}

double SteelFiberCompositeBar::posEnvlpRotlim(double strain) {
  double strainLimit = POS_INF_STRAIN;

  if (strain <= rot1) return POS_INF_STRAIN;
  if (strain > rot1 && strain <= rot2 && E2 < 0.0)
    strainLimit = rot1 - mom1 / E2;
  if (strain > rot2 && E3 < 0.0) strainLimit = rot2 - mom2 / E3;

  if (strainLimit == POS_INF_STRAIN)
    return POS_INF_STRAIN;
  else if (posEnvlpStress(strainLimit) > 0)
    return POS_INF_STRAIN;
  else
    return strainLimit;
}

double SteelFiberCompositeBar::negEnvlpRotlim(double strain) {
  double strainLimit = NEG_INF_STRAIN;

  if (strain >= rot1) return NEG_INF_STRAIN;
  if (strain < rot1 && strain >= rot2 && E2 < 0.0)
    strainLimit = rot1 - mom1 / E2;
  if (strain < rot2 && E3 < 0.0) strainLimit = rot2 - mom2 / E3;

  if (strainLimit == NEG_INF_STRAIN)
    return NEG_INF_STRAIN;
  else if (negEnvlpStress(strainLimit) < 0)
    return NEG_INF_STRAIN;
  else
    return strainLimit;
}
