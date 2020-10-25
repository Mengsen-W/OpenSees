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
                                                                        
// $Revision: 1.0 $
// $Date: 2012-09-17 10:51:44 $
// $Source: /usr/local/cvs/OpenSees/SRC/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h,v $
                                                                        
                                                                        
// File: ~/system_of_eqn/linearSOE/sparseGEN/PFEMSolver_Mumps.h
//
// Written: Minjie 
// Created: Sep 17 2012
//
=======

// $Revision$
// $Date$


// Written: Minjie 
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

#include <PFEMSolver_Mumps.h>
#include <PFEMLinSOE.h>
#include <iostream>
#include <cmath>
<<<<<<< HEAD
#include <Timer.h>
#include <mpi.h>

PFEMSolver_Mumps::PFEMSolver_Mumps(int r, int e, int h, int s)
    :PFEMSolver(), theSOE(0), sid(), pid(), myid(0),
     relax(r), err(e), host(h), sym(s)
{
=======
#include <vector>
#include <Timer.h>
#ifdef _MUMPS
#include <mpi.h>
#endif
#include <elementAPI.h>

#ifdef _AMGCL
#include <amgcl/adapter/crs_tuple.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/solver/fgmres.hpp>
#include <amgcl/solver/lgmres.hpp>
#include <amgcl/solver/gmres.hpp>
#include <amgcl/solver/cg.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/relaxation/ilu0.hpp>
#include <amgcl/amg.hpp>
#include <amgcl/profiler.hpp>
#include <amgcl/adapter/zero_copy.hpp>
#endif

void* OPS_PFEMSolver_Mumps()
{
    int numdata = 1;
    int relax = 20;
    int err = 0;
    int sym = 0;
    int add = 0;
    int print = 0;
    double ptol = 1e-4;
    double Bitol = 1e-16;
    int maxiter = 100;

    while (OPS_GetNumRemainingInputArgs() > 0) {
        const char* opt = OPS_GetString();
        if (strcmp(opt, "-relax") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &relax) < 0) {
                    opserr << "WARNING: failed to get relax\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-err") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &err) < 0) {
                    opserr << "WARNING: failed to get err\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-sym") == 0) {

            sym = 1;

        } else if (strcmp(opt, "-print") == 0) {

            print = 1;

        } else if (strcmp(opt, "-added-mass") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &add) < 0) {
                    opserr << "WARNING: failed to get add\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-ptol") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &ptol) < 0) {
                    opserr << "WARNING: failed to get ptol\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-Bitol") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetDoubleInput(&numdata, &Bitol) < 0) {
                    opserr << "WARNING: failed to get Bitol\n";
                    return 0;
                }
            }

        } else if (strcmp(opt, "-pmaxiter") == 0) {

            if (OPS_GetNumRemainingInputArgs() > 0) {
                if (OPS_GetIntInput(&numdata, &maxiter) < 0) {
                    opserr << "WARNING: failed to get err\n";
                    return 0;
                }
            }
        }
    }

    PFEMSolver_Mumps* theSolver = new PFEMSolver_Mumps(relax,err,add,sym,print,
                                                       ptol,maxiter,Bitol);
    return new PFEMLinSOE(*theSolver);
}

PFEMSolver_Mumps::PFEMSolver_Mumps(int r, int e, int a, int s, int p,
                                   double tol, int niter, double bitol)
        :PFEMSolver(), theSOE(0),
#ifdef _MUMPS
         sid(),
#endif
         relax(r), err(e), add(a), sym(s), print(p),
         ptol(tol), Bitol(bitol), pmaxiter(niter)
{
#ifdef _MUMPS
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    // mumps id
    sid.job = JOB_INIT;
    sid.par = 1;
    sid.sym = sym;
    sid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&sid);

    sid.irn = 0;
    sid.jcn = 0;
<<<<<<< HEAD

    pid.job = JOB_INIT;
    pid.par = 1;
    pid.sym = sym;
    pid.comm_fortran = USE_COMM_WORLD;
    dmumps_c(&pid);

    // get process id
    MPI_Comm_rank(MPI_COMM_WORLD, &myid);

=======
#endif
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
}

PFEMSolver_Mumps::~PFEMSolver_Mumps()
{
<<<<<<< HEAD
=======
#ifdef _MUMPS
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    sid.job = JOB_END;
    dmumps_c(&sid);

    if(sid.irn != 0) delete [] sid.irn;
    if(sid.jcn != 0) delete [] sid.jcn;
<<<<<<< HEAD

    pid.job = JOB_END;
    dmumps_c(&pid);
=======
#endif
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
}

int
PFEMSolver_Mumps::solve()
{
<<<<<<< HEAD
    // Timer timer;
    // timer.start();
=======
#ifdef _MUMPS
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    cs* M = theSOE->M;
    cs* Gft = theSOE->Gft;
    cs* Git = theSOE->Git;
    cs* L = theSOE->L;
<<<<<<< HEAD
    cs* Qt = theSOE->Qt;
    Vector& Mhat = theSOE->Mhat;
=======
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    Vector& Mf = theSOE->Mf;
    Vector& X = theSOE->X;
    Vector& B = theSOE->B;
    ID& dofType = theSOE->dofType;
    ID& dofID = theSOE->dofID;
<<<<<<< HEAD
=======
    int stage = theSOE->stage;
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0

    int Msize = M->n;
    int Isize = Git->n;
    int Ssize = Msize-Isize;
    int Fsize = Mf.Size();
    int Psize = L->n;
<<<<<<< HEAD
    int Pisize = Mhat.Size();
    int size = X.Size();

    // numeric LU factorization of M
    // call mumps for factorization
    MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Isize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Ssize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Fsize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Psize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Pisize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&size, 1, MPI_INT, 0, MPI_COMM_WORLD);
    
    if(Msize > 0) {
	sid.job = JOB_FACTORIZATION;
	dmumps_c(&sid);
    }
    
    if(sid.info[0] != 0) {
	opserr<<"WARNING: failed to factorize -- PFEMSolver_Mumps::solve\n";
	return -1;
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"factorization time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();
    
    // structure and interface predictor : deltaV1 = M^{-1} * rsi
    Vector deltaV1;
    if(Msize > 0) {
	if(myid == 0) {
	    deltaV1.resize(Msize);
	    deltaV1.Zero();
	
	    // rsi
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 2) {
		    deltaV1(rowid+Ssize) = B(i);   // rsi
		} else if(rowtype == 0) {
		    deltaV1(rowid) = B(i);         // rsi
		}
	    }
	    double* deltaV1_ptr = &deltaV1(0);
	    sid.rhs = deltaV1_ptr;
	    sid.nrhs = 1;
	}
	sid.job = JOB_SOLUTION;
	dmumps_c(&sid);
	if(sid.info[0] != 0) {
	    opserr<<"WARNING: failed to solve predictor -- PFEMSolver_Mumps::solve\n";
	    return -1;
	}
    }
    
    // fluid predictor: deltaVf1 = Mf^{-1} * rf
    Vector deltaVf1;
    if(Fsize > 0) {
	if(myid == 0) {
	    deltaVf1.resize(Fsize);
	    deltaVf1.Zero();
	    
	    // rf
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 1) {
		    if(Mf(rowid) == 0) {
			opserr<<"WANING: Zero Mf at location "<<rowid<<" ";
			opserr<<" - PFEMLinSOE::solve()\n";
			return -1;
		    }
		    deltaVf1(rowid) = B(i)/Mf(rowid);         // rf
		}
	    }
	}
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"dVs1 = "<<deltaV1.Norm()<<"\n";
    // 	opserr<<"dVf1 = "<<deltaVf1.Norm()<<"\n";
    // 	opserr<<"predictor time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();

    // Mi^{-1}, Msi^{-1}
    cs* invMi1 = cs_spalloc(Isize, Isize, 1, 1, 1);
    cs* invMsi1 = cs_spalloc(Ssize, Isize, 1, 1, 1);
    if (Isize > 0) {
	Vector eyes(Msize);
	double* eyes_ptr = &eyes(0);

	for (int j=0; j<Isize; j++) {
	    if (myid == 0) {
		// rhs
		eyes.Zero();
		eyes(j+Ssize) = 1.0;
		
		// M^{-1}*eyes
		sid.rhs = eyes_ptr;
		sid.nrhs = 1;		
	    }
	    sid.job = JOB_SOLUTION;
	    dmumps_c(&sid);
	    if(sid.info[0] != 0) {
		opserr<<"WARNING: failed to solve invMi -- PFEMSolver_Mumps::solve\n";
		return -1;
	    }
	    if (myid == 0) {
		// copy
		for (int i=0; i<Msize; i++) {
		    if (eyes(i) != 0.0) {
			if (i>=Ssize) {
			    cs_entry(invMi1,i-Ssize,j,eyes(i));
			} else {
			    cs_entry(invMsi1,i,j,eyes(i));
			}
		    }
		}

	    }
	}

    }
    cs* invMi = cs_compress(invMi1);
    cs* invMsi = cs_compress(invMsi1);
    cs_spfree(invMi1);
    cs_spfree(invMsi1);
    opserr<<"invMi = "<<cs_norm(invMi)<<"\n";

    // Gi, Mf^{-1}*Gf
    cs* Gi = 0;
    cs* Gf = 0;
    if(Fsize > 0) {
	if(myid == 0) {
	    Gi = cs_transpose(Git, 1);
	    Gf = cs_transpose(Gft, 1);
	    for(int j=0; j<Psize; j++) {
		for(int k=Gf->p[j]; k<Gf->p[j+1]; k++) {
		    int i = Gf->i[k];
		    double& x = Gf->x[k];
		    x /= Mf(i);
		}
	    }
	}
    }

    // solve for pressure
    Vector deltaP;
    if(Psize>0) {
	// assembled format
	ICNTL(pid,5,0);

	// input matrix is centralized on the host
	ICNTL(pid,18,0);

	// workspace relaxation: 20%
	ICNTL(pid,14,relax);

	// dense right hand side
	ICNTL(pid,20,0);

	// centralized right hand side
	ICNTL(pid,21,0);

	// error messages
	ICNTL(pid,1,err);
	ICNTL(pid,2,err);
	ICNTL(pid,3,err);
	ICNTL(pid,4,err);

	// get S
	cs* S = 0;
	if(myid == 0) {
	    deltaP.resize(Psize);
	    deltaP.Zero();
	    double* deltaP_ptr = &deltaP(0);

	    // Gft*deltaVf1
	    if(Fsize > 0) {
		double* deltaVf1_ptr = &deltaVf1(0);
		cs_gaxpy(Gft, deltaVf1_ptr, deltaP_ptr);
	    }

	    // Git*deltaVi1
	    if(Isize > 0) {
		double* deltaVi1_ptr = &deltaV1(0) + Ssize;
		cs_gaxpy(Git, deltaVi1_ptr, deltaP_ptr);
	    }

	    // rp-Git*deltaVi1-Gft*deltaVf1
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype == 3) {             // pressure
		    deltaP(rowid) = B(i)-deltaP(rowid);   // rp-Git*deltaVi1-Gft*deltaVf1
		}
	    }

	    // S = L + Git*Mi{-1}*Gi + Gft*Mf{-1}*Gf
	    if(Isize > 0) {
		cs* S1 = cs_multiply(Git, invMi);
		S = cs_multiply(S1, Gi);
		cs_spfree(S1);
	    }
	    if(Fsize > 0) {
		cs* S1 = cs_multiply(Gft, Gf);
		if(S == 0) {
		    S = S1;
		} else {
		    cs* S2 = cs_add(S, S1, 1.0, 1.0);
		    cs_spfree(S);
		    cs_spfree(S1);
		    S = S2;
		}
	    }

	    // S
	    if(S == 0) {
		S = L;
	    } else {
		cs* S1 = cs_add(S, L, 1.0, 1.0);
		cs_spfree(S);
		S = S1;
	    }

	    // set S
	    pid.n = S->n;
	    pid.nz = S->nzmax;
	    pid.a = S->x;
	    pid.irn = new int[pid.nz];
	    pid.jcn = new int[pid.nz];
	    for(int j=0; j<pid.n; j++) {
		for(int k=S->p[j]; k<S->p[j+1]; k++) {
		    pid.irn[k] = S->i[k]+1;
		    pid.jcn[k] = j+1;
		}
	    }

	    pid.rhs = deltaP_ptr;
	    pid.nrhs = 1;
	}
	// timer.pause();
	// if(myid == 0) {
	//     opserr<<"pressure setup time = "<<timer.getReal()<<"\n";
	// }
	
	// timer.start();
	
	// solve
	pid.job = JOB_SOLVE;
	dmumps_c(&pid);
	if(pid.info[0] != 0) {
	    opserr<<"WARNING: failed to solve pressure -- PFEMSolver_Mumps::solve\n";
	    return -1;
	}

	// release
	if (myid == 0) {
	    if(S != L) cs_spfree(S);
	    delete [] pid.irn;
	    delete [] pid.jcn;
	}
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"dP = "<<deltaP.Norm()<<"\n";
    // 	opserr<<"pressure  time = "<<timer.getReal()<<"\n";
    // }

    //timer.start();
    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    Vector deltaV;
    if(myid == 0) {
	deltaV.resize(Msize);
	deltaV.Zero();
    }
    if(Isize > 0) {
	if(myid == 0) {
	    // Gi*deltaP
	    Vector Gip(Isize);
	    double* Gip_ptr = &Gip(0);
	    if(Psize > 0) {
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Gi, deltaP_ptr, Gip_ptr);

		// Msi^{-1}*Gi*deltaP
		if(Ssize > 0) {
		    Vector vs(Ssize);
		    double* vs_ptr = &vs(0);
		    cs_gaxpy(invMsi, Gip_ptr, vs_ptr);
		    for(int i=0; i<Ssize; i++) {
			deltaV(i) += vs(i);
		    }
		}

		// Mi^{-1}*Gi*deltaP
		Vector vi(Isize);
		double* vi_ptr = &vi(0);
		cs_gaxpy(invMi, Gip_ptr, vi_ptr);
		for(int i=0; i<Isize; i++) {
		    deltaV(i+Ssize) = vi(i);
		}

	    }
	    cs_spfree(Gi);
	}
    }
    if(myid == 0) {
	deltaV += deltaV1;
    }
    cs_spfree(invMi);
    cs_spfree(invMsi);
    opserr<<"dVs = "<<deltaV.Norm()<<"\n";

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    Vector deltaVf(Fsize);
    if(Fsize > 0) {
        if(Psize > 0) {
	    if(myid == 0) {
		deltaVf.resize(Fsize);
		deltaVf.Zero();
		double* deltaVf_ptr = &deltaVf(0);
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Gf, deltaP_ptr, deltaVf_ptr);
		deltaVf += deltaVf1;
		cs_spfree(Gf);
	    }
        }
    }
    opserr<<"dVf = "<<deltaVf.Norm()<<"\n";

    // deltaPi = Mhatd^{-1}*(Qt*deltaP - rpi)
    Vector deltaPi;
    if(Pisize > 0) {
	if(myid == 0) {
	    deltaPi.resize(Pisize);
	    deltaPi.Zero();
	    
	    // Qt*deltaP
	    if(Psize > 0) {
		double* deltaPi_ptr = &deltaPi(0);
		double* deltaP_ptr = &deltaP(0);
		cs_gaxpy(Qt, deltaP_ptr, deltaPi_ptr);   // Qt*deltaP
	    }

	    // rpi
	    for(int i=0; i<size; i++) {        // row
		int rowtype = dofType(i);      // row type
		int rowid = dofID(i);          // row id
		if(rowtype != 4) continue;     // rpi
		if(Mhat(rowid) == 0.0) {
		    opserr<<"Zero Mhat at location "<<rowid<<" ";
		    opserr<<" - PFEMLinSOE::solve()\n";
		    return -1;
		}
		deltaPi(rowid) = (B(i)-deltaPi(rowid))/Mhat(rowid); // rpi
	    }
	}
        
    }
    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"corrector  time = "<<timer.getReal()<<"\n";
    // }
    // timer.start();
    // copy to X
    if(myid == 0) {
	X.Zero();

	for(int i=0; i<size; i++) {            // row
	    int rowtype = dofType(i);          // row type
	    int rowid = dofID(i);
	    if(rowtype == 0) {
		X(i) = deltaV(rowid);            
	    } else if(rowtype == 2) {
		X(i) = deltaV(rowid+Ssize);
	    } else if(rowtype == 1) {
		X(i) = deltaVf(rowid);
	    } else if(rowtype == 3) {
		X(i) = deltaP(rowid);
	    } else if(rowtype == 4) {
		X(i) = deltaPi(rowid);
	    }

	}
    }

=======
    int size = X.Size();

    // numeric LU factorization of M
    if(Msize > 0 && stage!=2) {
        sid.job = JOB_FACTORIZATION;
        dmumps_c(&sid);

        if(sid.info[0] != 0) {
            opserr<<"WARNING: failed to factorize -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }

    // structure and interface predictor : deltaV1 = M^{-1} * rsi
    std::vector<double> deltaV1;
    if(Msize > 0) {
        deltaV1.assign(Msize, 0.0);
    }
    if (Msize>0 && stage!=2) {
        // rsi
        for (int i = 0; i < size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if (rowtype == 2) {
                deltaV1[rowid + Ssize] = B(i);   // rsi
            } else if (rowtype == 0) {
                deltaV1[rowid] = B(i);         // rsi
            }
        }
        sid.rhs = &deltaV1[0];
        sid.nrhs = 1;
        sid.job = JOB_SOLUTION;
        dmumps_c(&sid);
        if (sid.info[0] != 0) {
            opserr << "WARNING: failed to solve predictor -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }

    // fluid predictor: deltaVf1 = Mf^{-1} * rf
    std::vector<double> deltaVf1;
    if (Fsize > 0) {
        deltaVf1.assign(Fsize,0.0);
    }
    if(Fsize>0 && stage!=2) {

        // rf
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 1) {
                if(Mf(rowid) == 0) {
                    opserr<<"WANING: Zero Mf at location "<<rowid<<" ";
                    opserr<<" - PFEMLinSOE::solve()\n";
                    return -1;
                }
                deltaVf1[rowid] = B(i)/Mf(rowid);         // rf
            }
        }
    }

    // Gi, Gf
    cs* Gi = 0;
    cs* Gf = 0;
    if(Fsize>0 && (stage==0||stage==2)) {
        Gf = cs_transpose(Gft, 1);
    }
    if(Isize>0 && (stage==0||stage==2)) {
        Gi = cs_transpose(Git, 1);
    }

    // solve for pressure
    std::vector<double> deltaP, rhsP;
    if(Psize>0) {
        deltaP.assign(Psize, 0.0);
        rhsP.assign(Psize, 0.0);
    }
    if(Psize>0 && (stage==2||stage==0)) {

        // S = L + Git*Mi{-1}*Gi + Gft*Mf{-1}*Gf
        cs* S = 0;

        // Gft*deltaVf1
        if(Fsize > 0 && stage==0) {
            cs_gaxpy(Gft, &deltaVf1[0], &rhsP[0]);
        }

        // Git*deltaVi1
        if(Isize > 0 && stage==0) {
            cs_gaxpy(Git, &deltaV1[0]+Ssize, &rhsP[0]);
        }

        // rp-Git*deltaVi1-Gft*deltaVf1
        for(int i=0; i<size; i++) {        // row
            int rowtype = dofType(i);      // row type
            int rowid = dofID(i);          // row id
            if(rowtype == 3) {             // pressure
                rhsP[rowid] = B(i)-rhsP[rowid];   // rp-Git*deltaVi1-Gft*deltaVf1
            }
        }

        // Git*Mi{-1}*Gi
        if (Isize > 0 && stage==0) {
            std::vector<int> irhs_ptr(Msize+1,1);
            std::vector<int> irhs_row(Isize*Isize);
            std::vector<double> rhs_val(Isize*Isize);

            for (int j=0; j<Isize; ++j) {
                for (int i=0; i<Isize; ++i) {
                    irhs_row[Isize*j+i] = i+Ssize+1;
                }
                irhs_ptr[j+1+Ssize] = (j+1)*Isize+1;
            }

            ICNTL(sid,20,3);
            ICNTL(sid,30,1);
            sid.nz_rhs = (int)rhs_val.size();
            sid.nrhs = Msize;
            sid.irhs_ptr = &irhs_ptr[0];
            sid.irhs_sparse = &irhs_row[0];
            sid.rhs_sparse = &rhs_val[0];

            sid.job = JOB_SOLUTION;
            dmumps_c(&sid);
            if(sid.info[0] != 0) {
                opserr<<"info[1] = "<<sid.info[0]<<", info[2] ="<<sid.info[1]<<"\n";
                opserr<<"WARNING: failed to solve Bi -- PFEMSolver_Mumps::solve\n";
                return -1;
            }
            ICNTL(sid,20,0);
            ICNTL(sid,30,0);

            for (int j=Isize; j<Msize; ++j) {
                irhs_ptr[j+1] -= 1;
            }
            for (int k=0; k<(int)irhs_row.size(); ++k) {
                irhs_row[k] -= 1;
                irhs_row[k] -= Ssize;
            }

            cs* Bi1 = cs_spalloc(Isize,Isize,1,1,1);
            int numignore = 0;
            for (int j=0; j<Isize; ++j) {
                for (int i=0; i<Isize; ++i) {
                    if (fabs(rhs_val[Isize*j+i]) > Bitol) {
                        cs_entry(Bi1, i, j,rhs_val[Isize*j+i]);
                    } else {
                        numignore++;
                    }
                }
            }
            cs* Bi = cs_compress(Bi1);
            cs_spfree(Bi1);

            cs* S1 = cs_multiply(Git, Bi);
            S = cs_multiply(S1, Gi);
            cs_spfree(S1);
            cs_spfree(Bi);

        } else if (Isize > 0 && stage==2) {

            // get Mi
            Vector Mi;
            Mi.resize(Isize);
            Mi.Zero();

            for (int i = 0; i < size; ++i) {
                int coltype = dofType(i);
                int colid = dofID(i);
                if (coltype == 2) {
                    int cid = colid + Ssize;
                    for (int k = M->p[cid]; k < M->p[cid+1]; ++k) {
                        Mi(colid) += M->x[k];
                    }
                }
            }

            // Gi*Mi{-1}*Gi
            for(int j=0; j<Psize; j++) {
                for(int k=Gi->p[j]; k<Gi->p[j+1]; k++) {
                    Gi->x[k] /= Mi(Gi->i[k]);
                }
            }
            cs* S1 = cs_multiply(Git, Gi);
            if(S == 0) {
                S = S1;
            } else {
                cs* S2 = cs_add(S, S1, 1.0, 1.0);
                cs_spfree(S);
                cs_spfree(S1);
                S = S2;
            }
        }

        // Gft*Mf{-1}*Gf
        if(Fsize > 0) {
            for(int j=0; j<Psize; j++) {
                for(int k=Gf->p[j]; k<Gf->p[j+1]; k++) {
                    Gf->x[k] /= Mf(Gf->i[k]);
                }
            }
            cs* S1 = cs_multiply(Gft, Gf);
            if(S == 0) {
                S = S1;
            } else {
                cs* S2 = cs_add(S, S1, 1.0, 1.0);
                cs_spfree(S);
                cs_spfree(S1);
                S = S2;
            }
        }

        // S
        if(S == 0) {
            S = L;
        } else {
            cs* S1 = cs_add(S, L, 1.0, 1.0);
            cs_spfree(S);
            S = S1;
        }

        // solve
        if (S->nzmax > 0) {
#ifdef _AMGCL
            // solve
            amgcl::profiler<> prof;
            typedef
            amgcl::make_solver<
                    amgcl::amg<
                            amgcl::backend::builtin<double>,
                            amgcl::coarsening::smoothed_aggregation,
                            amgcl::relaxation::spai0
                    >,
                    amgcl::solver::lgmres<amgcl::backend::builtin<double> >
            > Solver;


            // parameter
            Solver::params prm;
            prm.solver.tol = ptol;
            prm.solver.maxiter = pmaxiter;

            // setup
            prof.tic("setup");
            std::vector<std::ptrdiff_t> ptr(Psize+1), num(S->nzmax);
            for (int i=0; i<S->nzmax; ++i) {
                num[i] = S->i[i];
            }
            for (int i = 0; i < Psize+1; ++i) {
                ptr[i] = S->p[i];
            }
            double* val = &(S->x[0]);
            Solver solve(amgcl::adapter::zero_copy(Psize,&ptr[0],&num[0],val),prm);
            prof.toc("setup");

            // solve
            int iters;
            double error;
            prof.tic("solve");
            std::tie(iters, error) = solve(rhsP, deltaP);
            prof.toc("solve");

            if (print) {
                std::cout << solve << std::endl;
                std::cout << "iters: " << iters << std::endl
                          << "error: " << error << std::endl
                          << prof << std::endl;
            }

            if (iters>=pmaxiter && error>ptol) {
                opserr<<"WARNING: failed to solve pressure\n";
                return -1;
            }
#endif
        }

        // release
        if(S != L) cs_spfree(S);
    }

    // structure and interface corrector : deltaV = deltaV1 + M^{-1}*G*deltaP
    std::vector<double> deltaV;
    if (Msize > 0) {
        deltaV.assign(Msize, 0.0);
    }
    if(Isize > 0 && stage==0) {
        // Gi*deltaP
        if(Psize > 0) {
            cs_gaxpy(Gi, &deltaP[0], &deltaV[0]+Ssize);
        }

        sid.rhs = &deltaV[0];
        sid.nrhs = 1;
        sid.job = JOB_SOLUTION;
        dmumps_c(&sid);
        if(sid.info[0] != 0) {
            opserr<<"WARNING: failed to solve corrector -- PFEMSolver_Mumps::solve\n";
            return -1;
        }
    }
    for (int i = 0; i < Msize; ++i) {
        deltaV[i] += deltaV1[i];
    }

    // fluid corrector: deltaVf = deltaVf1 + Mf^{-1}*Gf*deltaP
    std::vector<double> deltaVf;
    if(Fsize > 0) {
        deltaVf.assign(Fsize, 0.0);
    }
    if(Fsize > 0 && stage==0) {
        if(Psize > 0) {
            cs_gaxpy(Gf, &deltaP[0], &deltaVf[0]);
        }
    }
    for (int i=0; i<Fsize; ++i) {
        deltaVf[i] = deltaVf[i] + deltaVf1[i];
    }

    // delete Gi, Gf
    if(Fsize>0 && (stage==0||stage==2)) {
        cs_spfree(Gf);
    }
    if(Isize>0 && (stage==0||stage==2)) {
        cs_spfree(Gi);
    }

    // copy to X
    X.Zero();

    for(int i=0; i<size; i++) {            // row
        int rowtype = dofType(i);          // row type
        int rowid = dofID(i);
        if(rowtype == 0) {
            X(i) = deltaV[rowid];
        } else if(rowtype == 2) {
            X(i) = deltaV[rowid+Ssize];
        } else if(rowtype == 1) {
            X(i) = deltaVf[rowid];
        } else if(rowtype == 3) {
            X(i) = deltaP[rowid];
        }
    }
#endif
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    return 0;
}

int PFEMSolver_Mumps::setSize()
{
<<<<<<< HEAD
    // Timer timer;
    // timer.start();

=======
#ifdef _MUMPS
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    // assembled format
    ICNTL(sid,5,0);

    // input matrix is centralized on the host
    ICNTL(sid,18,0);

    // workspace relaxation: 20%
    if(relax <= 0) {
<<<<<<< HEAD
	relax = 20;
=======
        relax = 20;
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    }
    ICNTL(sid,14,relax);

    // dense right hand side
    ICNTL(sid,20,0);

    // centralized right hand side
    ICNTL(sid,21,0);

    // No error messages
    if(err < 0) err = 0;
    ICNTL(sid,1,err);
    ICNTL(sid,2,err);
    ICNTL(sid,3,err);
    ICNTL(sid,4,err);

    // check M size
    cs* M = theSOE->M;
    int Msize = M->n;
<<<<<<< HEAD
    MPI_Bcast(&Msize, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(Msize <= 0) return 0;

    // host
    if(myid == 0) {
	sid.n = M->n;
	sid.nz = M->nzmax;
	sid.a = M->x;
	if(sid.irn != 0) delete [] sid.irn;
	if(sid.jcn != 0) delete [] sid.jcn;
	sid.irn = new int[sid.nz];
	sid.jcn = new int[sid.nz];

	for(int j=0; j<sid.n; j++) {
	    for(int k=M->p[j]; k<M->p[j+1]; k++) {
		sid.irn[k] = M->i[k]+1;
		sid.jcn[k] = j+1;
	    }
	}
=======
    if(Msize <= 0) return 0;

    // analysis
    sid.n = M->n;
    sid.nz = M->nzmax;
    sid.a = M->x;
    if(sid.irn != 0) delete [] sid.irn;
    if(sid.jcn != 0) delete [] sid.jcn;
    sid.irn = new int[sid.nz];
    sid.jcn = new int[sid.nz];

    for(int j=0; j<sid.n; j++) {
        for(int k=M->p[j]; k<M->p[j+1]; k++) {
            sid.irn[k] = M->i[k]+1;
            sid.jcn[k] = j+1;
        }
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    }

    // call mumps
    sid.job = JOB_ANALYSIS;
    dmumps_c(&sid);
    if(sid.info[0] != 0) {
<<<<<<< HEAD
	opserr<<"WARNING: failed to analyze -- PFEMSolver_Mumps::setSize\n";
	return -1;
    }

    // timer.pause();
    // if(myid == 0) {
    // 	opserr<<"analysis time = "<<timer.getReal()<<"\n";
    // }
=======
        opserr<<"WARNING: failed to analyze -- PFEMSolver_Mumps::setSize\n";
        return -1;
    }
#endif
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
    return 0;
}

int
PFEMSolver_Mumps::sendSelf(int cTag, Channel &theChannel)
{
    // nothing to do
    return 0;
}

int
PFEMSolver_Mumps::recvSelf(int ctag,
<<<<<<< HEAD
		  Channel &theChannel, 
		  FEM_ObjectBroker &theBroker)
=======
                           Channel &theChannel,
                           FEM_ObjectBroker &theBroker)
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
{
    // nothing to do
    return 0;
}


<<<<<<< HEAD
int 
=======
int
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
PFEMSolver_Mumps::setLinearSOE(PFEMLinSOE& theSOE)
{
    this->theSOE = &theSOE;
    return 0;
}

<<<<<<< HEAD
=======
#ifdef _MUMPS
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
void
PFEMSolver_Mumps::ICNTL(DMUMPS_STRUC_C& id, int I, int val)
{
    id.icntl[I-1] = val;
}
<<<<<<< HEAD


=======
#endif
>>>>>>> ad2965e00858958011abb8d72d2ec3efc732a9a0
