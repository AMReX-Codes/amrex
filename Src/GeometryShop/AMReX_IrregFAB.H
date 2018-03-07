#ifndef _IRREGFAB_H___
#define _IRREGFAB_H___

#include <cmath>
#include <cstdlib>
#include "AMReX_IntVectSet.H"
#include "AMReX_VolIndex.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_BaseFab.H"
#include "AMReX_BaseIVFAB.H"
#include "AMReX_EBGraph.H"

namespace amrex
{

///
/**
   This is a data holder which holds data over the irregular
   cells in a box.
*/
  class IrregFAB: public BaseIVFAB<Real>
  {
  public:
    ///
    /**
       Default constructor.  Constructs an uninitialized IVFAB.  You must
       subsequently call {\bf define} before using this IVFAB.
    */
    IrregFAB();

    ///
    /**
       Defining constructor.  Specifies the irregular domain as the
       irregular cells within the region
       and the number of data components per VoF.  The irregular domain
       must lie completely within the EBGraph.  The
       contents are uninitialized.  Calls full define function.
    */
    IrregFAB(const Box&        a_region,
             const EBGraph&    a_ebgraph,
             const int&        a_nvarin);


    ///
    IrregFAB(const IntVectSet& a_ivs,
             const EBGraph&    a_ebgraph,
             const int&        a_nvarin):
      BaseIVFAB<Real>(a_ivs, a_ebgraph, a_nvarin)
      {
        m_region = a_ivs.minBox();
      }


    ///
    ~IrregFAB();

    ///
    /**
       Full define function.  Specifies the irregular domain  as the
       irregular cells within the regionand the
       number of data components per VoF.  The irregular domain
       must lie completely within the EBGraph.   The contents
       are uninitialized.  If it has previously been defined, the old definition
       data is overwritten and lost.
    */
      void define(const Box&        a_region,
                  const EBGraph&    a_ebgraph,
                  const int&        a_nvarin);


    ///invalid but necessary for leveldata to compile
    IrregFAB(const Box& a_region, int a_nVar)
      {
        amrex::Error("invalid constructor called for baseebcellfab");
      }

    /// some stuff to do arithmetic
    class arithOp
    {
    public:
      virtual void func(Real& dstVal, const Real& srcVal) = 0;
    };

    /// operator for +=
    class additionOp: public arithOp
    {
      void func(Real& a_dstVal, const Real& a_srcVal)
        {
          a_dstVal += a_srcVal;
        }
    };

    /// operator for -=
    class subtractionOp: public arithOp
    {
      void func(Real& a_dstVal, const Real& a_srcVal)
        {
          a_dstVal -= a_srcVal;
        }
    };

    /// operator for *=
    class multiplicationOp: public arithOp
    {
      void func(Real& a_dstVal, const Real& a_srcVal)
        {
          a_dstVal *= a_srcVal;
        }
    };

    /// operator for /=
    class divisionOp: public arithOp
    {
      void func(Real& a_dstVal, const Real& a_srcVal)
        {
          a_dstVal /= a_srcVal;
        }
    };

    ///apply any operator pointwise for two fabs.
    IrregFAB& applyOp(const Real& a_srcVal,
                      int a_srccomp,
                      int a_destcomp,
                      int a_numcomp,
                      IrregFAB::arithOp& a_op);

    ///apply any operator pointwise with input real number
    IrregFAB& applyOp(const IrregFAB& a_srcVal,
                      const Box& a_region,
                      int a_srccomp,
                      int a_destcomp,
                      int a_numcomp,
                      IrregFAB::arithOp& a_op);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    IrregFAB& operator+=(const IrregFAB& a_ebfab);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    IrregFAB& operator+=(const Real& a_ebfab);


    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    IrregFAB& operator-=(const IrregFAB& a_ebfab);

    ///
    /**
     */
    IrregFAB& operator-=(const Real& a_ebfab);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    IrregFAB& operator*=(const IrregFAB& a_IrregFAB);

    ///
    /**
     */
    IrregFAB& operator*=(const Real& a_ebfab);

    ///
    /**
       done over intersection of src, dst
       Both fabs need the same ncomp
    */
    IrregFAB& operator/=(const IrregFAB& a_ebfab);

    ///
    /**
     */
    IrregFAB& operator/=(const Real& a_ebfab);


    ///more flexible versions of stuff

    ///
    IrregFAB& plus(const IrregFAB& a_src,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

    ///
    IrregFAB& plus(const Real& a_src,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

    ///
    IrregFAB& plus(const IrregFAB& a_src,
                   const Box& a_region,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

    ///
    IrregFAB& minus(const IrregFAB& a_src,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);

    ///
    IrregFAB& minus(const Real& a_src,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);

    ///
    IrregFAB& mult(const Real& a_src,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);
    ///
    IrregFAB& divide(const Real& a_src,
                     int a_srccomp,
                     int a_destcomp,
                     int a_numcomp);

    ///
    IrregFAB& minus(const IrregFAB& a_src,
                    const Box& a_region,
                    int a_srccomp,
                    int a_destcomp,
                    int a_numcomp);


    ///
    IrregFAB& mult(const IrregFAB& a_src,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

    ///
    IrregFAB& mult(const IrregFAB& a_src,
                   const Box& a_region,
                   int a_srccomp,
                   int a_destcomp,
                   int a_numcomp);

    ///
    IrregFAB& divide(const IrregFAB& a_src,
                     int a_srccomp,
                     int a_destcomp,
                     int a_numcomp);

    ///
    IrregFAB& divide(const IrregFAB& a_src,
                     const Box& a_region,
                     int a_srccomp,
                     int a_destcomp,
                     int a_numcomp);
  protected:
    Box m_region;

  private:

    //disallowed for all the usual reasons
    void operator= (const IrregFAB& a_input);
    IrregFAB(const IrregFAB& a_input);

  };
}
#endif
