
#include <AMReX_EBFabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_FabArray.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB2.H>

namespace amrex
{

EBFArrayBoxFactory::EBFArrayBoxFactory (const Geometry& a_geom,
                                        const BoxArray& a_ba,
                                        const DistributionMapping& a_dm,
                                        const Vector<int>& a_ngrow, EBSupport a_support)
    : m_support(a_support),
      m_geom(a_geom),
      m_ebdc(std::make_shared<EBDataCollection>(a_geom,a_ba,a_dm,a_ngrow,a_support))
{}

EBFArrayBoxFactory::EBFArrayBoxFactory (const EB2::Level& a_level,
                                        const Geometry& a_geom,
                                        const BoxArray& a_ba,
                                        const DistributionMapping& a_dm,
                                        const Vector<int>& a_ngrow, EBSupport a_support)
    : m_support(a_support),
      m_geom(a_geom),
      m_ebdc(std::make_shared<EBDataCollection>(a_level,a_geom,a_ba,a_dm,a_ngrow,a_support))
{}

FArrayBox*
EBFArrayBoxFactory::create (const Box& box, int ncomps,
                            const FabInfo& info, int box_index) const
{
    if (m_support == EBSupport::none)
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared);        
    }
    else
    {
        const EBCellFlagFab& ebcellflag = m_ebdc->getMultiEBCellFlagFab()[box_index];
        return new EBFArrayBox(ebcellflag, box, ncomps);
    }
}

EBFArrayBoxFactory*
EBFArrayBoxFactory::clone () const
{
    return new EBFArrayBoxFactory(*this);
}

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const Geometry& a_geom, const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support)
{
    std::unique_ptr<EBFArrayBoxFactory> r;
    if (EB2::use_eb2) {
        const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
        const EB2::Level& eb_level = index_space.getLevel(a_geom);
        r.reset(new EBFArrayBoxFactory(eb_level, a_geom, a_ba, a_dm, a_ngrow, a_support));
    } else {
        r.reset(new EBFArrayBoxFactory(a_geom, a_ba, a_dm, a_ngrow, a_support));
    }
    return r;
}

}
