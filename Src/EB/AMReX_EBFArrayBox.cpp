
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

EBFArrayBox::EBFArrayBox ()
    : FArrayBox(),
      m_ebisbox(),
      m_ebcellflag(nullptr)
{
}

EBFArrayBox::EBFArrayBox (const EBISBox& ebisBox, const EBCellFlagFab& ebflag,
                          const Box& bx, int ncomps)
    : FArrayBox(bx, ncomps),
      m_ebisbox(ebisBox),
      m_ebcellflag(&ebflag)
{
    const Box& sect = amrex::enclosedCells(bx) & ebisBox.getRegion();
    m_type = ebflag.getType(sect);
}

EBFArrayBox::~EBFArrayBox ()
{

}

const EBCellFlagFab&
getEBCellFlagFab (const FArrayBox& fab)
{
    const EBFArrayBox* ebfab = static_cast<EBFArrayBox const*>(&fab);
    BL_ASSERT(ebfab);
    return ebfab->getEBCellFlagFab();
}

const EBISBox&
getEBISBox (const FArrayBox& fab)
{
    const EBFArrayBox* ebfab = static_cast<EBFArrayBox const*>(&fab);
    BL_ASSERT(ebfab);
    return ebfab->getEBISBox();
}

}
