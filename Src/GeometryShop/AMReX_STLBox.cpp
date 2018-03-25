#include "AMReX_STLBox.H"

namespace amrex
{


  STLBox::STLBox(shared_ptr<STLMesh>    a_stlmesh,
                 const Box&             a_region,
                 const Box          &   a_domain,
                 const RealVect&        a_origin,
                 const RealVect&        a_dx)
  {
    SetMeshBox(a_stlmesh,a_region,a_domain,a_origin,a_dx); // just set some stuff in other method
  }

  void STLBox::SetMeshBox(shared_ptr<STLMesh>   a_stlmesh,
                          const Box&             a_region,
                          const Box          &   a_domain,
                          const RealVect&        a_origin,
                          const RealVect&        a_dx)
  {
    // set data, used in construction or to reset data later (e.g. after default construction)
    m_msh     = a_stlmesh;
    m_region  = a_region;
    m_domain  = a_domain;
    m_origin  = a_origin;
    m_dx      = a_dx;
  }

}
