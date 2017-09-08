
#include <AMReX_EBFaceFlag.H>

namespace amrex {

constexpr std::array<std::array<int,3>,3> EBFaceFlag::pos_ngbr;

EBFaceFlagFab::EBFaceFlagFab (const Box& bx, int, bool, bool)
    : m_box(bx),
      m_fabs {{ AMREX_D_DECL(BaseFab<EBFaceFlag>(surroundingNodes(grow(bx,-1),0)),
                             BaseFab<EBFaceFlag>(surroundingNodes(grow(bx,-1),1)),
                             BaseFab<EBFaceFlag>(surroundingNodes(grow(bx,-1),2))) }}
{
}

EBFaceFlagFab::~EBFaceFlagFab ()
{
}

}
