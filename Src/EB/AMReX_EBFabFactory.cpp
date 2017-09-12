
#include <AMReX_EBFabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>

namespace amrex
{

EBFArrayBoxFactory::EBFArrayBoxFactory (const EBLevel& a_eblevel)
    : m_eblevel(a_eblevel)
{}

FArrayBox*
EBFArrayBoxFactory::create (const Box& box, int ncomps,
                            const FabInfo& info, int box_index) const
{
    const auto& ebisl = this->getEBISL();
    const auto& eblevel = this->getEBLevel();
    if (ebisl.isDefined())
    {
        const EBISBox& ebisBox = ebisl[box_index];
        const EBCellFlagFab& ebcellflag = eblevel.getMultiEBCellFlagFab()[box_index];
        return new EBFArrayBox(ebisBox, ebcellflag, box, ncomps);
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
