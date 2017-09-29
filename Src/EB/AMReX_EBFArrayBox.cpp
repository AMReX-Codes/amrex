
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex {

EBFArrayBox::EBFArrayBox ()
    : FArrayBox(),
      m_ebcellflag(nullptr)
{
}

EBFArrayBox::EBFArrayBox (const EBCellFlagFab& ebflag, const Box& bx, int ncomps)
    : FArrayBox(bx, ncomps),
      m_ebcellflag(&ebflag)
{
    BL_ASSERT(ebflag.box().strictly_contains(amrex::enclosedCells(bx)));
    const Box& sect = amrex::enclosedCells(bx) & ebflag.getRegion();
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

}
