#ifdef CH_LANG_CC
/*
 *      _______              __
 *     / ___/ /  ___  __ _  / /  ___
 *    / /__/ _ \/ _ \/  V \/ _ \/ _ \
 *    \___/_//_/\___/_/_/_/_.__/\___/
 *    Please refer to Copyright.txt, in Chombo's root directory.
 */
#endif

#if defined(CH_Darwin) && defined(__GNUC__) && ( __GNUC__ == 3 )
// deal with the broken isnan()/isinf() in GCC on MacOS
#include <unistd.h>
#define _GLIBCPP_USE_C99 1
#endif

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
using namespace std;

#include "GeometryService.H"
#include "GeometryShop.H"
#include "Moments.H"
#include "ConstrainedLS.H"
#include "PolyGeom.H"
#include "RealVect.H"
#include "parstream.H"

#include "NamespaceHeader.H"

ConstrainedLS::ConstrainedLS()
  :
  m_nbound(0),
  m_residual(LARGEREALVAL)
{
}

ConstrainedLS::LSResult ConstrainedLS::qrSolution(Real ** a,
                                                  Vector<Real> & x,
                                                  Vector<Real> & rhs,
                                                  Real         & resq)
{
  const int nCol = x.size();
  const int nRow = rhs.size();

  resq=-2.0;               // underdetermined
  if (nRow < nCol) return UNDERDETERMINED;
  resq=-1.0;               // singular

  //Rotates  into upper triangular form.
  for (int j = 0; j<nCol; ++j)
    {
      //   Find constants for rotation and diagonal entry.
      Real sq=0.0;
      for (int i = j; i<nRow; ++i)
        {
          sq+= a[i][j]*a[i][j];
        }
      if (sq == 0.0)
      {
        resq = -1.0; //singular
        return SINGULAR;
      }
      Real qv1= a[j][j] >= 0 ? -sqrt(sq) : sqrt(sq);
      Real u1=a[j][j] - qv1;
      a[j][j]=qv1;
      int j1=j + 1;
      //  Rotate remaining columns of sub-matrix.
      for (int jj = j1; jj<nCol ; ++jj)
        {
          Real dot=u1*a[j][jj];
          for (int i = j1; i<nRow ; ++i)
            {
              dot+=a[i][jj]*a[i][j];
            }
          Real constant=dot/Abs(qv1*u1);
          for (int i = j1; i<nRow ; ++i)
            {
              a[i][jj]-= constant*a[i][j]; // todo pointer
            }
          a[j][jj] -= constant*u1;
        }

      //  Rotate  rhs  vector.
      Real dot=u1*rhs[j];
      for (int i=j1; i<nRow; ++i)
        {
          dot+=rhs[i]*a[i][j];
        }

      Real constant=dot/Abs(qv1*u1);
      rhs[j]-=constant*u1;
      for (int i=j1; i<nRow; ++i)
        {
          rhs[i] -= constant*a[i][j];
        }
    }

  //  Solve triangular system by back-substitution.
  for (int ii = 0; ii<nCol; ++ii)
    {
      int i=nCol-(ii+1);
      Real sum=rhs[i];
      for (int j = i+1; j<nCol; ++j)
        {
          sum=sum - a[i][j]*x[j];
        }
      if (a[i][i] == 0.0) return SINGULAR;
      x[i]=sum/a[i][i];
    }
  // Find residual in overdetermined case.
  resq=0.0;
  for (int i = nCol; i<nRow; ++i)
    {
      resq+=rhs[i]*rhs[i];
    }
  return SUCCESS;
}

Real ConstrainedLS::getResidual() const
{
  return m_residual;
}


int ConstrainedLS::numberActiveConstraints() const
{
  return m_nbound;
}

Vector<ConstrainedLS::Bound> ConstrainedLS::getConstraints() const
{
  return m_boundState;
}


ConstrainedLS::LSResult ConstrainedLS::solveUnconstrained(Vector<Real>       & a_x,
                                           Real**               a_A,
                                           const Vector<Real> & a_rhs)
{
  Vector<Real>  lowBound(a_x.size());
  lowBound.assign(-HUGE);
  Vector<Real> highBound(a_x.size());
  highBound.assign(HUGE);
  return solveBoundConstrained(a_x,a_A,a_rhs,lowBound,highBound);
}

bool ConstrainedLS::boundsConsistent(const Vector<Real> & a_loBound,
                                     const Vector<Real> & a_hiBound) const
{
  bool retval = true;
  Real maxDiff = -LARGEREALVAL;
  for (int jbound = 0; jbound < a_loBound.size(); ++jbound)
  {
    retval &= (a_loBound[jbound] <= a_hiBound[jbound]);
    maxDiff = Max(maxDiff, (a_hiBound[jbound] - a_loBound[jbound]));
  }
  retval &= (maxDiff >= 0.0);
  return retval;
}

ConstrainedLS::LSResult ConstrainedLS::solveBoundConstrained(Vector<Real>      & x,
                                                             Real**              A,
                                                             const Vector<Real>& rhs,
                                                             const Vector<Real>& lowerBound,
                                                             const Vector<Real>& upperBound)
{
  const int NO_IFROM5 = -1234567;
  LSResult lastResultQR;
  // Note: in the past this eps actually mattered a lot in domains with flat sections
  // it got changed from 1.e-16
  const Real eps = 1.e-13;
  int numColsA = x.size();
  int numRowsA = rhs.size();
  int n = numColsA;
  int m = numRowsA;
  Real** act = NULL;
  allocArray(numRowsA,numColsA,act);
  Vector<Real> actResid(numRowsA);

  Vector<int> istate(x.size()+1);

  Real sj = 0.0; //todo: init correctly

  /*  Step 1.  Initialize everything--
      active and bound sets, initial values, etc.*/
  int mm = m <= n ? m : n;  // min m,n
  int jj = -1;       // index of last successful entry to active set
  int ifrom5 =NO_IFROM5;   // indicator step 6 entered from step 5 (sign
                    // gives bound, value gives index)
  bool from5Low = false;
  int iact = -1; //todo: check this
  int key = 1;   //todo: hardwire
  Vector<Real> w(n);
  Vector<Real> zz(n);
  vector<bool> isBoundLow(n);


  /*  Check bounds*/
  CH_assert(lowerBound.size() == upperBound.size());
  if (!boundsConsistent(lowerBound,upperBound))
    {
      pout() << "Inconsistent bounds in BVLS constrained least squares algorithm"<<endl;;
      return INCONSISTENT_BOUNDS;
    }

  m_nbound = 0;
  int nact = n-m_nbound;


  /* This initialization is only the warm start with all variables
     assumed active */

  for (int k = 0; k < numColsA ; ++k)
    {
      if (lowerBound[k] < 0.0 && upperBound[k] > 0.0)
        {
          x[k] = 0.0;
        }
      else if (lowerBound[k] == -HUGE)
        {
          x[k]=upperBound[k] - eps;
        }
      else if (upperBound[k] == HUGE)
        {
          x[k] = lowerBound[k] + eps;
        }
      else
        {
          x[k]=(lowerBound[k]+upperBound[k])/2.0;
        }
      istate[k] = k;
      isBoundLow[k]=false;
    }
  istate[numColsA] = 0; // m_nbound = 0


  // Compute bnorm, the norm of the data vector b, for reference.
  // todo:  Lapack?
  Real bsq=0.0;
  for (int i = 0; i < numRowsA ; ++ i)
    {
      bsq+= rhs[i]*rhs[i];
    }
  Real bnorm=sqrt(bsq);
  for (int iLoopA = 1; iLoopA < 3*n; ++iLoopA)
    {
      // Step 2: Initialize the negative gradient vector w
      Real obj = 0.0;
      for (int j=0; j < n; ++j)
        {
          w[j]=0.0;
        }
      for (int i=0; i<m; ++i)
        {
          Real ri=rhs[i];
          for (int j=0; j<n; ++j)
            {
              ri -= A[i][j]*x[j];
            }
          obj+= ri*ri;   // obj = || a.x - b ||.
          for (int j=0; j<n; ++j)
            {
              w[j] += A[i][j]*ri;
            }
          actResid[i]=ri;  //The residual vector is stored in the
                          //mm+1'st column of act(*,*).
                          // todo: seems kludgy
        }
      // Converged?  Stop if the misfit << || b ||,
      // or if all components are active (unless this is the
      // first iteration from a warm start).

      if ((sqrt(obj) <= bnorm*eps) || (iLoopA > 1 && m_nbound == 0))
        {
          istate[numColsA]=m_nbound;
          w[0]=sqrt(obj);
          m_residual=sqrt(obj);
          freeArray(numRowsA,numColsA,act);
          return lastResultQR;
        }

      // Add the contribution of the active components back into the residual.
        for (int k = m_nbound; k < n ; ++k)
          {
            int j=istate[k];
            for (int i=0; i<m; ++i)
              {
                actResid[i] += A[i][j]*x[j];
              }

          }

        //The first iteration in a warm start requires immediate qr.

        if (!(iLoopA == 1 && key != 0))
          {
            //Find the bound element that most wants to be active.
            findBound: Real worst=0.0;
            int it=1;
            for (int j=0; j<m_nbound; ++j)
              {
                int ks=istate[j];
                Real bad=isBoundLow[ks] ? -w[ks] : w[ks];
                if (bad < worst)
                  {
                    it=j;
                    worst=bad;
                    iact=ks; //todo: scope
                  }
              }

            // Test whether the Kuhn-Tucker condition is met.

            if (worst >= 0.0 )
              {
                istate[n]=m_nbound;
                w[0]=sqrt(obj);
                m_residual=sqrt(obj);
                freeArray(numRowsA,numColsA,act);
                return SUCCESS;
              }


         // The component  x(iact)  is the one that most wants to become active.
         // If the last successful change in the active set was to move x(iact)
         // to a bound, don't let x(iact) in now: set the derivative of the
         // misfit with respect to x(iact) to zero and return to the Kuhn-Tucker
         // test.
         if ( iact == jj )
           {
             w[jj]=0.0;
             goto findBound; //todo: came from old Fortran
           }

         // Step 5. Undo the effect of the new (potentially)
         // active variable on the residual vector.
         Real bound = isBoundLow[iact] ? lowerBound[iact] : upperBound[iact];
         for (int i=0; i<m; ++i)
           {
             actResid[i]+= bound*A[i][iact];
           }

         //  Set flag ifrom5, indicating that Step 6 was entered from Step 5.
         //  This forms the defined but not usedasis of a test for instability: the gradient
         //  calculation shows that x(iact) wants to join the active
         //  set; if
         //  qr puts x(iact) beyond the bound from which it came, the gradient
         //   calculation was in error and the variable should not have been
         //   introduced.

         ifrom5=istate[it]; //todo?
         from5Low=isBoundLow[ifrom5];

         //  Swap the indices (in istate) of the new active variable and the
         //   rightmost bound variable; `unbind' that location by decrementing
         //   m_nbound.

         istate[it]=istate[m_nbound-1]; //todo?
         m_nbound--;
         nact++;
         istate[m_nbound]=iact; //todo?
         isBoundLow[istate[m_nbound]] = false;

         if (mm < nact)
           {
             pout() << "Too many free variables in BVLS constrained least squares algorithm!";
             return UNDERDETERMINED;
           }
         }

         //Step 6.
         //Load array  act  with the appropriate columns of  a  for qr.  For
         //added stability, reverse the column ordering so that the most
         //recent addition to the active set is in the last column.  Also
         //copy the residual vector from act(., mm1) into act(.,
         //mm1+1).
         bool doQR = true;
         while (doQR)
           {
             Vector<Real> actResidCopy(numRowsA);
             zz.resize(nact);
             actResidCopy = actResid;
             // prepare for qr using active problem
             for (int i=0; i<m; ++i)
               {
                 for (int k=m_nbound; k<n; ++k)
                   {
                     int j=istate[k];
                     act[i][nact-1-k+m_nbound]=A[i][j]; //todo?
                   }
               }
             Real resq;
             lastResultQR = qrSolution(act,zz,actResidCopy,resq);

             //  Test for linear dependence in qr, and for an instability that moves
             //   the variable just introduced away from the feasible region
             //   (rather than into the region or all the way through it).
             //   In either case, remove the latest vector introduced from the
             //   active set and adjust the residual vector accordingly.
             //   Set the gradient component (w(iact)) to zero and return to
             //   the Kuhn-Tucker test.

             if (ifrom5 !=NO_IFROM5)
             {
               if (resq < 0.0
                   || (!from5Low && (zz[nact-1] > upperBound[iact]))
                   || ( from5Low && (zz[nact-1] < lowerBound[iact])))
                 {
                   m_nbound++;
                   isBoundLow[istate[m_nbound]] = (x[iact]-upperBound[iact] < 0.);

                   nact--;
                   for (int i=0; i<m; ++i)
                     {
                       actResid[i] -= x[iact]*A[i][iact]; //todo? iact
                     }
                   ifrom5=NO_IFROM5;
                   w[iact]=0.0;
                   goto findBound;   // came from old-school Fortran
                 }
               /*  If Step 6 was entered from Step 5 and we are here, a new variable
                   has been successfully introduced into the active set; the last
                   variable that was fixed at a bound is again permitted to become
                   active.*/
               jj=-1;
             }
             ifrom5=NO_IFROM5;
             int k1;
             bool foundInfeasible = false;
             // Step 7.  Check for strict feasibility of the new qr solution.
             for (int k=0; k<nact; ++k)
               {
                 k1=k;
                 int j=istate[k+m_nbound];
                   if (zz[nact-1-k] < lowerBound[j] ||
                       zz[nact-1-k] > upperBound[j])
                     {
                       foundInfeasible = true;
                       break;
                     }
               }
             if (! foundInfeasible)
               {
                 for (int k=0; k<nact; ++k)
                   {
                     int j=istate[k+m_nbound];
                     x[j]=zz[nact-1-k];
                   }
                 break; // get out of doQR to top of main loop
               }

             //  Steps 8, 9.
             Real alpha=2.0;
             Real alf=alpha;
             for (int k=k1; k<nact; ++k)//todo:
               {
                 int j=istate[k+m_nbound];
                 if (zz[nact-1-k] > upperBound[j])
                   {
                     alf=(upperBound[j]-x[j])/(zz[nact-1-k]-x[j]);
                   }
                 if (zz[nact-1-k] < lowerBound[j])
                   {
                     alf=(lowerBound[j]-x[j])/(zz[nact-1-k]-x[j]);
                   }
                 if (alf < alpha)
                   {
                     alpha=alf;
                     jj=j;
                     sj=(zz[nact-1-k]-lowerBound[j]) >= 0.0 ? 1. : -1.;
                   }
               }

             //  Step 10
             for (int k=0; k<nact; ++k)
               {
                 int j=istate[k+m_nbound];
                 x[j] += alpha*(zz[nact-1-k]-x[j]);
               }

             /*  Step 11.
                 Move the variable that determined alpha to the appropriate bound.
                 (jj is its index; sj is + if zz(jj)> upperBound(jj)] - if zz(jj)<bl(jj) ).
                 If any other component of  x  is infeasible at this stage, it must
                 be due to roundoff.  Bind every infeasible component and every
                 component at a bound to the appropriate bound.  Correct the
                 residual vector for any variables moved to bounds.  Since at least
                 one variable is removed from the active set in this step, Loop B
                 (Steps 6-11) terminates after at most  nact  steps.*/

             int noldb=m_nbound;
             for (int k=0; k<nact; ++k)
               {
                 int j=istate[k+noldb];
                 if (((upperBound[j]-x[j]) <= 0.0) ||
                     ((j==jj) && (sj > 0.0)))
                   {
                     x[j]=upperBound[j];
                     istate[k+noldb]=istate[m_nbound];
                     istate[m_nbound]=j;
                     m_nbound++;
                     for (int i=0; i<m; ++i)
                       {
                         actResid[i]-= upperBound[j]*A[i][j];
                       }
                   }
                 else if (((x[j]-lowerBound[j]) <= 0.0) ||
                          ((j == jj) && (sj < 0.0)))
                   {
                     x[j]=lowerBound[j];
                     istate[k+noldb]=istate[m_nbound];
                     istate[m_nbound]=j;
                     isBoundLow[j]=true;
                     m_nbound++;
                     for (int i=0; i<m ; ++i)
                       {
                         actResid[i]-= lowerBound[j]*A[i][j];
                       }
                   }
               }
             nact=n - m_nbound;
             // If there are still active variables left repeat the qr;
             doQR = (nact > 0 );
           } //doQR
    } //LoopA
  pout() << "BVLS constrained least squares algorithm failed to converge" <<endl;
  freeArray(numRowsA,numColsA,act);  //
                                     //todo: this isn't done anywhere
                                     //where it will be called
  return UNCONVERGED;
}



  void ConstrainedLS::allocArray(const int& rows,
                                 const int& cols,
                                 Real**& A)
{
  A = new Real* [rows];

  for (int i = 0; i < rows;i++)
    {
      A[i] = new Real [cols];
      Real* scanA = A[i];

      for (int j = 0; j < cols; j++)
        {
          *(scanA++) = 0.0;
        }
    }
}

void ConstrainedLS::freeArray(const int& rows,
                              const int& cols,
                              Real**& A)
{
  for (int i = 0; i < rows; i++)
    {
      delete[] A[i];
    }

  delete[] A;
}



#include "NamespaceFooter.H"
