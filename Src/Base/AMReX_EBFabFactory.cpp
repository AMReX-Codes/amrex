
#include <AMReX_FabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
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
        const EBISBox& ebisBox = m_ebisl[box_index];
        return new EBFArrayBox(ebisBox, box, ncomps);
    }
    else
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared);        
    }
}

}
