#ifndef _BASEEBBC_H_
#define _BASEEBBC_H_

#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_VolIndex.H"
#include "AMReX_EBISLayout.H"
#include "AMReX_EBCellFAB.H"
#include "AMReX_EBFluxFAB.H"
#include "AMReX_Stencils.H"
#include "AMReX_EBLevelGrid.H"
#include "AMReX_BaseBCFuncEval.H"
#include "AMReX_LayoutData.H"


namespace amrex
{
///
/**
 */
  class ConductivityBaseEBBC
  {
  public:
    ///
    /**
     */
    ConductivityBaseEBBC()
      {
        m_coefSet  = false;
        m_onlyHomogeneous = true;
      }

    virtual ~ConductivityBaseEBBC()
      {
      }


    /// the operator calls this
    /**
     */
    virtual void setCoef(EBLevelGrid                                &  a_eblg,
                         Real                                       &  a_beta,
                         Real                                       &  a_dx,
                         const RealVect                             &  a_probLo,
                         shared_ptr<FabArray<EBFluxFAB> >           &  a_bcoe)
      {
        m_coefSet = true;
        m_beta    = a_beta;
        m_eblg    = a_eblg;
        m_bcoe    = a_bcoe;
        m_dx      = a_dx;
        m_probLo  = a_probLo;
      }


    ///
    /**
       Return a pointer to the homogenous flux stencil for the boundary condition.
       contribution.   In the case where there in no contribution in the homogeneous case
       (ie. homogeneous Neumann) return NULL.
    */
    virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)= 0;

    ///
    /**
       add change in lphi due to inhomogeneous eb flux
    */
    virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                             const EBCellFAB&              a_phi,
                             const vector<VolIndex>&       a_vofsToChange,
                             const MFIter  &               a_mfi,
                             const Real&                   a_factor,
                             const bool&                   a_useHomogeneous) = 0;

    ///
    void
    setValue(Real a_value)
      {
        m_onlyHomogeneous = false;
        m_isFunction = false;
        m_value = a_value;
        m_func = shared_ptr<BaseBCFuncEval>();
      }
    ///
    void
    setFunction(shared_ptr<BaseBCFuncEval> a_func)
      {
        m_value = 12345.6789;
        m_func = a_func;

        m_onlyHomogeneous = false;
        m_isFunction = true;
      }

    Real bcvaluefunc(const RealVect      & a_point)
      {
        Real retval;
        if(m_onlyHomogeneous)
        {
          retval = 0;
        }
        else if(m_isFunction)
        {
          retval = m_func->value(a_point, 0);
        }
        else
        {
          retval = m_value;
        }
        return retval;
      }
  protected:

    EBLevelGrid                              m_eblg;
    bool                                     m_coefSet;
    Real                                     m_beta;
    Real                                     m_dx;
    RealVect                                 m_probLo;
    shared_ptr<FabArray<EBFluxFAB> >         m_bcoe;
    Real m_value;
    shared_ptr<BaseBCFuncEval>        m_func;
    bool m_onlyHomogeneous;
    bool m_isFunction;

  };
///
/**
 */
  class ConductivityBaseEBBCFactory
  {
  public:
    ConductivityBaseEBBCFactory()
      {
        m_onlyHomogeneous = true;
      }

    virtual ~ConductivityBaseEBBCFactory()
      {
      }

    ///
    /**
     */
    virtual ConductivityBaseEBBC* new_object_ptr() = 0;

    ///
    void
    setValue(Real a_value)
      {
        m_onlyHomogeneous = false;
        m_isFunction = false;
        m_value = a_value;
        m_func = shared_ptr<BaseBCFuncEval>();
      }
    ///
    void
    setFunction(shared_ptr<BaseBCFuncEval> a_func)
      {
        m_onlyHomogeneous = false;
        m_value = 12345.6789;
        m_func = a_func;
        m_isFunction = true;
      }

    ///
    ConductivityBaseEBBC* create()
      {
        ConductivityBaseEBBC* newBC = new_object_ptr();
        if (!m_onlyHomogeneous)
        {
          if (m_isFunction)
          {
            newBC->setFunction(m_func);
          }
          else
          {
            newBC->setValue(m_value);
          }
        }
        else
        {
          newBC->setValue(0.0);
        }
        return static_cast<ConductivityBaseEBBC*>(newBC);
      }
  protected:

    bool m_onlyHomogeneous;
    bool m_isFunction;
    Real m_value;
    shared_ptr<BaseBCFuncEval>        m_func;

  };
}

#endif
