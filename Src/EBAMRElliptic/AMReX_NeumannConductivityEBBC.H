#ifndef _NEUMANNCONDUCTIVITYEBBC_H_
#define _NEUMANNCONDUCTIVITYEBBC_H_

#include "AMReX_ConductivityBaseEBBC.H"
#include "AMReX_EBCellFAB.H"

namespace amrex
{
  ///
  /**
   */
  class NeumannConductivityEBBC: public ConductivityBaseEBBC
  {
  public:
    ///
    virtual void applyEBFlux(EBCellFAB&                    a_lphi,
                             const EBCellFAB&              a_phi,
                             const vector<VolIndex>&       a_vofsToChange,
                             const MFIter  &               a_mfi,
                             const Real&                   a_factor,
                             const bool&                   a_useHomogeneous);

    ///
    NeumannConductivityEBBC():ConductivityBaseEBBC()
      {
      }
    
    virtual LayoutData<BaseIVFAB<VoFStencil> >* getFluxStencil(int ivar)
      {
        return NULL;
      }
    
    ///
    virtual ~NeumannConductivityEBBC()
      {
      }
    
    
  private:

  };
    
  ///
  class NeumannConductivityEBBCFactory: public ConductivityBaseEBBCFactory
  {
  public:
    ///
    NeumannConductivityEBBCFactory()
      {
      }
    
    ///
    virtual ~NeumannConductivityEBBCFactory()
      {
      }
    
    
    ///
    virtual ConductivityBaseEBBC* new_object_ptr()
      {
        NeumannConductivityEBBC* fresh = new NeumannConductivityEBBC();

        return static_cast<ConductivityBaseEBBC*>(fresh);
      }

    
  private:

    
  };
}    
#endif
