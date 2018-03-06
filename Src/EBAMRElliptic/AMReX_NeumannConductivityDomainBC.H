#ifndef _NEUMANNCONDUCTIVITYDOMAINBC_H___
#define _NEUMANNCONDUCTIVITYDOMAINBC_H___

#include "AMReX_ConductivityBaseDomainBC.H"

namespace amrex
{
///
/**
   Class to do Dirichlet boundary conditions on the domain for the conductivity operator.
*/
  class NeumannConductivityDomainBC: public ConductivityBaseDomainBC
  {
  public:
    NeumannConductivityDomainBC();

    virtual ~NeumannConductivityDomainBC();

    ///
    virtual void 
    fillPhiGhost(FArrayBox&     a_state,
                 const Box&     a_valid,
                 bool           a_homogeneous);


    ///this is for when the EB and the domain intersect
    virtual void getFaceFlux(Real&                 a_faceFlux,
                             const VolIndex&       a_vof,
                             const MFIter    &     a_mfi,
                             const EBCellFAB&      a_phi,
                             const int&            a_idir,
                             const Side::LoHiSide& a_side,
                             const bool&           a_useHomogeneous);

    ///
    void
    getFaceGradPhi(Real&                 a_faceGradPhi,
                   const FaceIndex&      a_face,
                   const MFIter    &     a_mfi,
                   const EBCellFAB&      a_phi,
                   const int&            a_idir,
                   const Side::LoHiSide& a_side,
                   const bool&           a_useHomogeneous);

  private:


  };

  class NeumannConductivityDomainBCFactory: public ConductivityBaseDomainBCFactory
  {
  public:
    NeumannConductivityDomainBCFactory()
      {
      }

    virtual ~NeumannConductivityDomainBCFactory()
      {
      }

    ConductivityBaseDomainBC* new_object_ptr()
      {
        NeumannConductivityDomainBC* newBC = new NeumannConductivityDomainBC();
        return static_cast<ConductivityBaseDomainBC*>(newBC);
      }

  };
}
#endif
