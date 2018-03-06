// We can assume that template class T has null construction.

#ifndef _LINEARSOLVER_H___
#define _LINEARSOLVER_H___

#include "AMReX_REAL.H"
#include "AMReX_Box.H"
#include <cmath>

namespace amrex
{
///
/**
   Operator class for Linear solver for solving L(phi) = rhs.
*/
  template <class T>
  class LinearOp
  {
  public:
    ///
    virtual ~LinearOp()
      {
      }

    ///
    /**
       Say you are  solving L(phi) = rhs.   Make a_lhs = L(a_phi) - a_rhs.   If a_homogeneous is true,
       evaluate the operator using homogeneous boundary conditions.
    */
    virtual void residual(  T& a_lhs, const T& a_phi, const T& a_rhs, bool a_homogeneous = false) = 0;

    ///
    /**
       Given the current state of the residual the correction, apply your preconditioner to a_cor.
    */
    virtual void preCond(   T& a_cor, const T& a_residual)       = 0;

    ///
    /**
       In the context of solving L(phi) = rhs, set a_lhs = L(a_phi).  If a_homogeneous is true,
       evaluate the operator using homogeneous boundary conditions.
    */
    virtual void applyOp(   T& a_lhs, const T& a_phi, bool a_homogeneous = false) = 0;

    ///
    /**
       Creat data holder a_lhs that mirrors a_rhs.   You do not need to copy the data of a_rhs,
       just  make a holder the same size.
    */
    virtual void create(    T& a_lhs, const T& a_rhs) = 0;

    ///
    /**
       Opposite of create -- perform any operations required before lhs goes
       out of scope. In general, the only time this needs to be defined in
       a derived class is if the create() function called new. Otherwise, the
       default no-op function is sufficient.
    */
    virtual void clear(T& a_lhs)
      {
      }

    ///
    /**
       Set a_lhs  equal to a_rhs.
    */
    virtual void assign(    T& a_lhs, const T& a_rhs)       = 0;

    ///
    /**
       Compute and return the dot product of a_1 and a_2.   In most contexts, this
       means return the sum over all data points of a_1*a_2.
    */
    virtual Real dotProduct(const T& a_1, const T& a_2)     = 0;

    ///
    /**
       Increment by scaled amount (a_lhs += a_scale*a_x).
    */
    virtual void incr  (    T& a_lhs, const T& a_x, Real a_scale) = 0;

    ///
    /**
       Set input to a scaled sum (a_lhs = a_a*a_x + a_b*a_y).
    */
    virtual void axby(      T& a_lhs, const T& a_x, const T& a_y, Real a_a, Real a_b) = 0;

    ///
    /**
       Multiply the input by a given scale (a_lhs *= a_scale).
    */
    virtual void scale(     T& a_lhs, const Real& a_scale)  = 0;

    ///
    /**
       Return the norm of  a_rhs.
       a_ord == 0  max norm, a_ord == 1 sum(abs(a_rhs)), else, L(a_ord) norm.
    */
    virtual Real norm(const T& a_rhs, int a_ord) = 0;

    ///
    /**
       Set a_lhs to zero.
    */
    virtual void setToZero(T& a_lhs) = 0;


  };

  ///
  /**
   Generic linear solver template.   BiCGStab and others are built on top of this.
  */
  template <class T>
  class LinearSolver
  {
  public:
    virtual ~LinearSolver()
      {
      }

    ///
    /**
       reset whether the solver is homogeneous.
    */
    virtual void setHomogeneous(bool a_homogeneous) = 0;

    ///
    /**
       Define the operator and whether it is a homogeneous solver or not.
       The LinearSolver does not take over ownership of this a_operator object.
       It does not call delete on it when the LinearSolver is deleted.  It is
       meant to  be like a late-binding reference. If you created a_operator
       with new, you should call delete on it after LinearSolver is deleted
       if you want to avoid memory leaks.

    */
    virtual void define(LinearOp<T>* a_operator, bool a_homogeneous = false) = 0;

    ///
    /**
       Solve L(phi) = rhs (phi = L^-1 (rhs)).
    */
    virtual void solve(T& a_phi, const T& a_rhs) = 0;

    /// Set a convergence metric, along with solver tolerance, if desired
    /**
       Default implementation does nothing, since there are probably
       cases (liked direct solves), where this has no real meaning.
    */
    virtual void setConvergenceMetrics(Real a_metric, Real a_tolerance)
      {
      }
  };

}

#endif /*_LINEARSOLVER_H_*/
