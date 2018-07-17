#ifndef AMREX_EBCELLFAB_H_
#define AMREX_EBCELLFAB_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_BaseEBCellFAB.H"
#include "AMReX_FArrayBox.H"
#include "AMReX_EBISBox.H"
#include "AMReX_Box.H"

namespace amrex
{
///
/**
   EBCellFAB is a BaseEBCellFAB<Real>
   which includes arithmetic functions.
*/
  class EBCellFAB: public BaseEBCellFAB<Real>
  {
  public:
    typedef Real value_type;

    ///
    EBCellFAB();

    ///
    /**
       Box going into this needs to be cell-centered.
       Data will exist over the surrounding nodes of the box.
    */
    EBCellFAB(const EBISBox& a_ebisBox,
              const Box& a_region, int a_nVar);

  

    ///
    void define(const EBISBox& a_ebisBox,
                const Box& a_region, int a_nVar);


    ///
    virtual ~EBCellFAB();

    ///
    const FArrayBox& getFArrayBox() const;

 
    ///
    /**
       Does a dynamic cast  of the BaseFab<Real>
       to a  FArrayBox for arithmetic purposes.
    */
    FArrayBox& getFArrayBox();

    ///
    BaseIVFAB<Real>& getMultiValuedFAB()
      {
        return m_irrFAB;
      }

    ///
    const BaseIVFAB<Real>& getMultiValuedFAB() const
      {
        return m_irrFAB;
      }
    ///
    /**
       Negate the values in all locations
    */
    EBCellFAB& negate(void);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    EBCellFAB& operator+=(const EBCellFAB& a_ebfab);


    EBCellFAB&
    axby(const EBCellFAB& a_X, 
         const EBCellFAB& a_Y,
         const Real& a_A, const Real& a_B);

    ///
    /**
       done over intersection of src, dst
    */
    EBCellFAB& plus(const EBCellFAB& a_ebfab,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp,
                    Real scale = 1.0);
    ///
    /**
       done over a_region
    */
    EBCellFAB& plus(const EBCellFAB& a_ebfab,
                    const Box& a_region,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp,
                    Real scale = 1.0);


    //for fabarray
    EBCellFAB& plus(const EBCellFAB& a_ebfab,
                    const Box& a_region1,
                    const Box& a_region2,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp)
      {
        Box region = a_region1 & a_region2;
        this->plus(a_ebfab, region, a_srccomp, a_destcomp, a_numcomp);
        return *this;
      }
    


    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    EBCellFAB& operator-=(const EBCellFAB& a_ebfab);

    ///
    /**
       done over intersection of src, dst
    */
    EBCellFAB& minus(const EBCellFAB& a_ebfab,
                     int a_srccomp,
                     int a_destcomp,
                     int a_numcomp);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    EBCellFAB& operator*=(const EBCellFAB& a_ebfab);

    ///
    /**
       done over intersection of src, dst
    */
    EBCellFAB&  mult(const EBCellFAB& a_src,
                     int a_srccomp,
                     int a_destcomp,
                     int a_numcomp);
    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    EBCellFAB& operator/=(const EBCellFAB& a_ebfab);


    ///
    /**
       done over intersection of src, dst
    */
    EBCellFAB& divide(const EBCellFAB& a_ebfab,
                      int a_srccomp,
                      int a_destcomp,
                      int a_numcomp);

    ///
    /**
     */
    EBCellFAB& operator+=(const Real& a_scalar);

    ///
    /**
     */
    EBCellFAB& operator-=(const Real& a_scalar);

    ///
    /**
     */
    EBCellFAB& operator*=(const Real& a_scalar);

    ///
    /**
     */
    EBCellFAB& mult(Real a_scalar);

    ///
    /**
     */
    EBCellFAB& operator/=(const Real& a_scalar);

    /// multiply each cell value by its volume fraction
    void kappaWeight();

    ///
    /**
       Current FAB = a_X FAB * a_A  +  a_Y FAB * a_B.
    */
    // EBCellFAB& axby(const EBCellFAB& a_X, const EBCellFAB& a_Y,
    //                 const Real& a_A, const Real& a_B,
    //                 const int& a_destComp,const int& a_xComp,const int& a_yComp);

    //! Returns the maximum value in the FAB for the given component.
    Real max(int a_comp = 0) const;

    //! Returns the minimum value in the FAB for the given component.
    Real min(int a_comp = 0) const;


    //needed for leveldata to compile
    //but very invalid
    EBCellFAB(const Box& a_region, int a_nComp)
      {
        amrex::Error("ebcelfab 1 invalid operator");
      }

    /// good for making temps.
    void clone(const EBCellFAB& a_arg);

  private:
    //disallowed for all the usual reasons
    EBCellFAB(const EBCellFAB& ebcin);
    void operator=(const EBCellFAB& fabin);

  };


}
#endif
