#ifndef AMREX_EBFACEFAB_H_
#define AMREX_EBFACEFAB_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_BaseEBFaceFAB.H"
#include "AMReX_REAL.H"
#include "AMReX_FArrayBox.H"

namespace amrex
{

  ///
  /**
     EBFaceFAB is a BaseEBFaceFAB<Real>
     which includes arithmetic functions.
  */
  class EBFaceFAB: public BaseEBFaceFAB<Real>
  {
  public:
    ///
    EBFaceFAB();
       
    ///
    /**
       Box going into this needs to be cell-centered.
       Data will exist over the surrounding nodes of the box.
    */
    EBFaceFAB(const EBISBox& a_ebisBox,
              const Box& a_region,
              int a_iDir, int a_nVar);
       
       
    void define(const EBISBox&  a_ebisBox,
                const Box& a_region,
                int a_iDir, int a_nComp);
    ///
    virtual ~EBFaceFAB();
       
    ///
    const FArrayBox& getFArrayBox() const;
       
    ///
    /**
       Does a dynamic cast of the BaseFab<Real>
       to a  FArrayBox for arithmetic purposes.
    */
    FArrayBox& getFArrayBox();
       
    ///
    /**
       Both fabs need the same ebisBox and region and number of components.
    */
    EBFaceFAB& operator+=(const EBFaceFAB& a_ebfab);
       
       
    ///
    /**
       Both fabs need the same ebisBox and region.
       negative component means nComp()
    */
    EBFaceFAB& plus(const EBFaceFAB& a_ebfab,
                    int a_srccomp = 0,
                    int a_destcomp= 0,
                    int a_numcomp = -1,
                    Real scale = 1.0);
       
       
    ///
    /**
       Both fabs need the same ebisBox and region and number of components.
    */
    EBFaceFAB& operator-=(const EBFaceFAB& a_ebfab);
       
    ///
    /**
       Both fabs need the same ebisBox and region.
       negative component means nComp()
    */
    EBFaceFAB& minus(const EBFaceFAB& a_ebfab,
                     int a_srccomp = 0,
                     int a_destcomp = 0,
                     int a_numcomp = -1,
                     Real scale = 1.0);
       
       
    ///
    /**
       Both fabs need the same ebisBox and region and number of components.
    */
    EBFaceFAB& operator*=(const EBFaceFAB& a_ebfab);
       
    ///
    /**
       Both fabs need the same ebisBox and region.
       negative component means nComp()
    */
    EBFaceFAB& mult(const EBFaceFAB& a_ebfab,
                    int a_srccomp = 0,
                    int a_destcomp = 0,
                    int a_numcomp = -1);
       
       
    ///
    /**
       Both fabs need the same ebisBox and region and number of components.
    */
    EBFaceFAB& operator/=(const EBFaceFAB& a_ebfab);
       
       
    ///
    /**
       Both fabs need the same ebisBox and region.
       negative component means nComp()
    */
    EBFaceFAB& divide(const EBFaceFAB& a_ebfab,
                      int a_srccomp = 0,
                      int a_destcomp = 0,
                      int a_numcomp= -1);
       
       
    ///
    /**
     */
    EBFaceFAB& operator+=(const Real& a_scalar);
       
    ///
    /**
     */
    EBFaceFAB& operator*=(const Real& a_scalar);
       
       
    //! Returns the maximum value for the given component in the FAB.
    Real max(int a_comp = 0) const;
       
    //! Returns the minimum value for the given component in the FAB.
    Real min(int a_comp = 0) const;
       
    //needed for leveldata to compile but invalid
    EBFaceFAB(const Box& a_region, int a_nComp)
      {
        amrex::Error("bogus constructor");
      }
  private:
    //verboten for all the usual reasons
    EBFaceFAB(const EBFaceFAB& ebcin);
    void operator=(const EBFaceFAB& ebcin);
       
  };
}
       
#endif
