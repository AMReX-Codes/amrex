#include "AMReX_EBFluxFactory.H"

namespace amrex
{
/***************/
  EBFluxFactory::~EBFluxFactory()
  {
  }
/***************/
  EBFluxFactory::EBFluxFactory(const EBISLayout& a_ebisl)
  {
    m_ebisl = a_ebisl;
  }
/***************/
  EBFluxFAB*
  EBFluxFactory::
  create(const Box& a_box, int a_ncomps, const FabInfo& info, int a_box_index) const
  {
    return new EBFluxFAB(m_ebisl[a_box_index], a_box, a_ncomps);
  }
/***************/
}
