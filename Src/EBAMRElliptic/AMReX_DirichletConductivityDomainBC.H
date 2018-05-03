#ifndef _DIRICHLETCONDUCTIVITYDOMAINBC_H___
#define _DIRICHLETCONDUCTIVITYDOMAINBC_H___

#include "AMReX_ConductivityBaseDomainBC.H"

namespace amrex
{

///
/**
   Class to do Dirichlet boundary conditions on the domain for the conductivity operator.
*/
  class DirichletConductivityDomainBC: public ConductivityBaseDomainBC
  {

  public:
    DirichletConductivityDomainBC();


    virtual ~DirichletConductivityDomainBC();


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
                   const FaceIndex&       a_face,
                   const MFIter    &     a_mfi,
                   const EBCellFAB&      a_phi,
                   const int&            a_idir,
                   const Side::LoHiSide& a_side,
                   const bool&           a_useHomogeneous);

  private:
  };

  ///factory to create the bcs.   After it is created you need to call fillCoefs.
  class DirichletConductivityDomainBCFactory: public ConductivityBaseDomainBCFactory
  {
  public:
    DirichletConductivityDomainBCFactory()
      {
      }

    virtual ~DirichletConductivityDomainBCFactory()
      {
      }

    ConductivityBaseDomainBC*
    new_object_ptr()
      {
        DirichletConductivityDomainBC* newBC = new DirichletConductivityDomainBC();
        return static_cast<ConductivityBaseDomainBC*>(newBC);
      }

  private:

  };

}
#endif
