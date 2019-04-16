
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
    BL_ASSERT(ebflag.box().contains(amrex::enclosedCells(bx)));
    const Box& sect = amrex::enclosedCells(bx) & ebflag.getRegion();
    m_type = ebflag.getType(sect);
}

EBFArrayBox::EBFArrayBox (EBFArrayBox const& rhs, MakeType make_type, int scomp, int ncomp)
    : FArrayBox(rhs, make_type, scomp, ncomp)
{
    m_type = rhs.m_type; // xxxxx TODO gpu
    m_ebcellflag = rhs.m_ebcellflag;
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
