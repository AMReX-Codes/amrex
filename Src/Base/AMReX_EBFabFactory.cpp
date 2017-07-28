
#include <AMReX_FabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBISLayout.H>

namespace amrex
{

FArrayBoxFactory::FArrayBoxFactory (const EBISLayout& a_ebisl)
    : m_ebisl(a_ebisl)
{}

FArrayBox*
FArrayBoxFactory::create (const Box& box, int ncomps,
                          const FabInfo& info, int box_index) const
{
    if (m_ebisl.isDefined())
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared);
    }
    else
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared);        
    }
}

}
