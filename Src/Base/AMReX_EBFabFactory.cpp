
#include <AMReX_EBFabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>

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
    if (m_ebisl.isDefined())
    {
        const EBISBox& ebisBox = ebisl[box_index];
        return new EBFArrayBox(ebisBox, box, ncomps);
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
