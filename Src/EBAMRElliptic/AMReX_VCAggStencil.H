#ifndef _VCAGGSTENCIL_H___
#define _VCAGGSTENCIL_H___

#include "AMReX_Stencils.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_REAL.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_AggStencil.H"


namespace amrex
{
  /// 
  /**
     specialization of aggstencil for 
     variable coefficient elliptic operators.
     Real* dataPtr(int dataType, int ivar);
  
     sten_t classes need the following functions
     srcIndex_t index(int isten)
     Real       weight(int isten)
  */
  class VCAggStencil: public AggStencil<EBCellFAB, EBCellFAB>
  {
  public:
    ///
    /**
       a_dstIndex is where the answers will go.
       a_srcStencil is the stencil into the src data
       a_phiData is to provide offsets for source data (can be a dummy but must be the right size)
       a_rhsData is to provide offsets for destin data (can be a dummy but must be the right size)
       ncomp is the number of components in the cache
       rhsData and phiData can have the wrong number of comps
    */
    VCAggStencil(const vector<shared_ptr<BaseIndex  > >   & a_dstVoFs,
                 const vector<shared_ptr<BaseStencil> >   & a_stencil,
                 const EBCellFAB                          & a_phiData,
                 const EBCellFAB                          & a_rhsData,
                 const EBCellFAB                          & a_relCoef,
                 const BaseIVFAB<Real>                    & a_alphaDiagWeight,
                 const int                                & a_ncomp);

    ///
    virtual ~VCAggStencil()
      {
      }

    ///
    virtual void   cachePhi(const EBCellFAB& a_phi) const;

    ///
    virtual void uncachePhi(      EBCellFAB& a_phi) const;

    ///
    /**
       for the given color
       phi := phi + relcoef*(lphi - rhs)
       where lphi = alpha*alphadiagweight*phi + beta*divF
       where divF is given by the input stencil above
    */
    virtual void relax(EBCellFAB              & a_phi,
                       const EBCellFAB        & a_rhs,
                       const EBCellFAB        & a_relCoef,
                       const BaseIVFAB<Real>  & a_alphaDiagWeight,
                       const Real             & a_alpha,
                       const Real             & a_beta,
                       const int              & a_varDest,
                       const IntVect          & a_color);

    ///
    /**
       lphi = alpha*alphadiagweight*phi + beta*divF
       where divF is given by the input stencil above
    */
    ///
    virtual void apply(EBCellFAB              & a_lphi,
                       const EBCellFAB        & a_phi,
                       const BaseIVFAB<Real>  & a_alphaDiagWeight,
                       const Real             & a_alpha,
                       const Real             & a_beta,
                       const int              & a_varDest,
                       const bool             & a_incrmentOnly);

  protected:
    vector<IntVect>                m_iv;
    int m_destVar;
    vector<access_t>               m_phiAccess;
    vector<access_t>               m_relAccess;

    vector<access_t>               m_alpAccess;

    vector<IntVect>                m_destIV;
    mutable vector< vector<Real> > m_cachePhi;

  private:
    /// disallowed operators.   Without code because Jeffster says that is better.

    ///weak construction is bad.
    VCAggStencil();

    ///deep copy for this object would kill performance
    VCAggStencil& operator=(const VCAggStencil& stenin);

    ///no reason to allow this one either.
    VCAggStencil(const VCAggStencil& stenin);

  };

}
#endif
