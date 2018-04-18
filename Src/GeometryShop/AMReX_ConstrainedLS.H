#ifndef AMREX_CONSTRAINEDLS_H_
#define AMREX_CONSTRAINEDLS_H_

#include "AMReX_Vector.H"
#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"


using amrex::Vector;
using amrex::Real;

class ConstrainedLS
{
public:
  enum Bound
  {
    LOWER_BOUND,
    UPPER_BOUND,
    UNCONSTRAINED
  };

  enum LSResult
  {
    SUCCESS,
    INCONSISTENT_BOUNDS,
    UNDERDETERMINED,
    SINGULAR,
    UNCONVERGED
  };

  ConstrainedLS();

  ///
  /** Solve the bound-constrained least squares
      Uses the Lawson-Hanson active set method
      as given by P.B. Stark and R.L. Parker,
      Bounded-Variable Least-Squares: An Algorithm
      and Applications, Computational Statistics, 10:129-141, 1995.
      Uses solveUnconstrained for sub-problems
  */
  LSResult solveBoundConstrained(Vector<Real>      & a_x,
                                 Real**              a_A,
                                 const Vector<Real>& a_rhs,
                                 const Vector<Real>& a_lowerBound,
                                 const Vector<Real>& a_upperBound);

  /** Solve the unconstrained least squares problem
   */
  LSResult solveUnconstrained(Vector<Real>       & a_x,
                              Real**             a_A,
                              const Vector<Real> & a_rhs);

  bool boundsConsistent(const Vector<Real> & a_lowBound,
                        const Vector<Real> & a_hiBound) const;

  /** Get a vector of Bound indicating constraint activity.
   */
  Vector<Bound> getConstraints() const;

  /** Query the number of active constraints.
      Bails if a LS problem has not been solved yet
  */
  int numberActiveConstraints() const;

  /** Get the LS residual
      Bails if a LS problem has not been solved yet
  */
  Real getResidual() const;

  /** Solve the least squares problem.
      Uses successive Householder rotations
  */
  LSResult qrSolution(Real **      a_A,
                      Vector<Real> & a_x,
                      Vector<Real> & a_b,
                      Real         & resq);




protected:
  Vector<Bound> m_boundState;
  int m_nbound;



private:
  Real m_residual;
  void allocArray(const int& rows,
                  const int& cols,
                  Real**& A);

  void freeArray(const int& rows,
                 const int& cols,
                 Real**& A);


};


#endif
