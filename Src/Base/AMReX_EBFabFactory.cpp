
#include <AMReX_EBFabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex
{

EBFArrayBoxFactory::EBFArrayBoxFactory (const EBLevel& a_eblevel)
    : m_eblevel(a_eblevel),
      m_ebisl(m_eblevel.getEBISL())
{}

FArrayBox*
EBFArrayBoxFactory::create (const Box& box, int ncomps,
                            const FabInfo& info, int box_index) const
{
    const auto& ebisl = this->getEBISLayout();
    const auto& eblevel = this->getEBLevel();
    if (m_ebisl.isDefined())
    {
        const EBISBox& ebisBox = ebisl[box_index];
        const EBFlagFab& ebflag = eblevel.Flags()[box_index];
        return new EBFArrayBox(ebisBox, ebflag, box, ncomps);
    }
    else
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared);        
    }
}

EBFArrayBoxFactory*
EBFArrayBoxFactory::clone () const
{
    return new EBFArrayBoxFactory(*this);
}

}
