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

// $Revision$
// $Date$
// $URL$

#ifndef HHTHSIncrReduct_TP_h
#define HHTHSIncrReduct_TP_h

// Written: Andreas Schellenberg (andreas.schellenberg@gmail.com)
// Created: 10/05
// Revision: A
//
// Description: This file contains the class definition for HHTHSIncrReduct_TP.
// HHTHSIncrReduct_TP is an algorithmic class for performing a transient analysis
// using the HHT integration scheme based on the trapezoidal rule. This is a special
// integrator for hybrid simulation where a reduction factor, that modifies the response
// increments, is introduced to enhance the smoothness of the convergence path and to
// reduce the possibility of spurious loading/unloading cycles during iteration.

#include <TransientIntegrator.h>

class DOF_Group;
class FE_Element;
class Vector;

class HHTHSIncrReduct_TP : public TransientIntegrator
{
public:
    // constructors
    HHTHSIncrReduct_TP();
    HHTHSIncrReduct_TP(double rhoInf, double reduct);
    HHTHSIncrReduct_TP(double alphaI, double alphaF,
        double beta, double gamma, double reduct);
    
    // destructor
    ~HHTHSIncrReduct_TP();
    
    // method to set up the system of equations
    int formUnbalance(void);
    
    // methods which define what the FE_Element and DOF_Groups add
    // to the system of equation object.
    int formEleTangent(FE_Element *theEle);
    int formNodTangent(DOF_Group *theDof);
    int formEleResidual(FE_Element *theEle);
    int formNodUnbalance(DOF_Group *theDof);
    
    int domainChanged(void);
    int newStep(double deltaT);
    int revertToLastStep(void);
    int update(const Vector &deltaU);
    int commit(void);
<<<<<<< HEAD
=======

    const Vector &getVel(void);
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    
    virtual int sendSelf(int commitTag, Channel &theChannel);
    virtual int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
    
    void Print(OPS_Stream &s, int flag = 0);
    
protected:
    
private:
    double alphaI;
    double alphaF;
    double beta;
    double gamma;
    double reduct;   // disp increment reduction factor 0<reduct<=1
    double deltaT;
    
    double c1, c2, c3;                      // some constants we need to keep
    double alphaM, alphaD, alphaR, alphaP;  // weighting factors we need to keep
    Vector *Ut, *Utdot, *Utdotdot;          // response quantities at time t
    Vector *U, *Udot, *Udotdot;             // response quantities at time t + deltaT
    Vector *scaledDeltaU;                   // scaled displacement increment
    Vector *Put;                            // unbalance at time t
};

#endif
