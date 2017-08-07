
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

EBFArrayBox::EBFArrayBox ()
    : FArrayBox(),
      m_ebisbox(),
      m_ebflag(nullptr)
{
}

EBFArrayBox::EBFArrayBox (const EBISBox& ebisBox, const EBFlagFab& ebflag,
                          const Box& box, int ncomps)
    : FArrayBox(box, ncomps),
      m_ebisbox(ebisBox),
      m_ebflag(&ebflag)
{
    const Box& sect = amrex::enclosedCells(box) & ebisBox.getRegion();
    m_type = ebflag.getType(sect);
}

EBFArrayBox::~EBFArrayBox ()
{

}

const EBFlagFab&
getEBFlagFab (const FArrayBox& fab)
{
    const EBFArrayBox* ebfab = static_cast<EBFArrayBox const*>(&fab);
    BL_ASSERT(ebfab);
    return ebfab->getEBFlagFab();
}

const EBISBox&
getEBISBox (const FArrayBox& fab)
{
    const EBFArrayBox* ebfab = static_cast<EBFArrayBox const*>(&fab);
    BL_ASSERT(ebfab);
    return ebfab->getEBISBox();
}

}
