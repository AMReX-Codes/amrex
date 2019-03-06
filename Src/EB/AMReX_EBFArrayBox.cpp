
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

#ifdef AMREX_USE_GPU
EBFArrayBox::EBFArrayBox (EBFArrayBox const& rhs, MakeType make_type)
    : FArrayBox(rhs, make_type)
{
    m_type = rhs.m_type; // xxxxx TODO gpu
    m_ebcellflag = rhs.m_ebcellflag;
}
#endif

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
