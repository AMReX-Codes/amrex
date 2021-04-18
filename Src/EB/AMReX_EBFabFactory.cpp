
#include <AMReX_EBFabFactory.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_EBFArrayBox.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_FabArray.H>

#include <AMReX_EB2_Level.H>
#include <AMReX_EB2.H>

namespace amrex
{

EBFArrayBoxFactory::EBFArrayBoxFactory (const EB2::Level& a_level,
                                        const Geometry& a_geom,
                                        const BoxArray& a_ba,
                                        const DistributionMapping& a_dm,
                                        const Vector<int>& a_ngrow, EBSupport a_support)
    : m_support(a_support),
      m_geom(a_geom),
      m_ebdc(std::make_shared<EBDataCollection>(a_level,a_geom,a_ba,a_dm,a_ngrow,a_support)),
      m_parent(&a_level)
{}

AMREX_NODISCARD
FArrayBox*
EBFArrayBoxFactory::create (const Box& box, int ncomps,
                            const FabInfo& info, int box_index) const
{
    if (m_support == EBSupport::none)
    {
        return new FArrayBox(box, ncomps, info.alloc, info.shared, info.arena);
    }
    else
    {
        const EBCellFlagFab& ebcellflag = m_ebdc->getMultiEBCellFlagFab()[box_index];
        return new EBFArrayBox(ebcellflag, box, ncomps, info.arena);
    }
}

AMREX_NODISCARD
FArrayBox*
EBFArrayBoxFactory::create_alias (FArrayBox const& rhs, int scomp, int ncomp) const
{
    if (m_support == EBSupport::none)
    {
        return new FArrayBox(rhs, amrex::make_alias, scomp, ncomp);
    }
    else
    {
        EBFArrayBox const& ebrhs = static_cast<EBFArrayBox const&>(rhs);
        return new EBFArrayBox(ebrhs, amrex::make_alias, scomp, ncomp);
    }
}

void
EBFArrayBoxFactory::destroy (FArrayBox* fab) const
{
    if (m_support == EBSupport::none)
    {
        delete fab;
    }
    else
    {
        EBFArrayBox* p = static_cast<EBFArrayBox*>(fab);
        delete p;
    }
}

AMREX_NODISCARD
EBFArrayBoxFactory*
EBFArrayBoxFactory::clone () const
{
    return new EBFArrayBoxFactory(*this);
}

bool
EBFArrayBoxFactory::isAllRegular () const noexcept
{
    return m_parent->isAllRegular();
}

EB2::IndexSpace const*
EBFArrayBoxFactory::getEBIndexSpace () const noexcept
{
    return (m_parent) ? m_parent->getEBIndexSpace() : nullptr;
}

int
EBFArrayBoxFactory::maxCoarseningLevel () const noexcept
{
    if (m_parent) {
        EB2::IndexSpace const* ebis = m_parent->getEBIndexSpace();
        return EB2::maxCoarseningLevel(ebis, m_geom);
    } else {
        return EB2::maxCoarseningLevel(m_geom);
    }
}

const DistributionMapping&
EBFArrayBoxFactory::DistributionMap () const noexcept
{
    return m_ebdc->getMultiEBCellFlagFab().DistributionMap();
}

const BoxArray&
EBFArrayBoxFactory::boxArray () const noexcept
{
    return m_ebdc->getMultiEBCellFlagFab().boxArray();
}

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support)
{
    const EB2::IndexSpace& index_space = EB2::IndexSpace::top();
    const EB2::Level& eb_level = index_space.getLevel(a_geom);
    return std::make_unique<EBFArrayBoxFactory>(eb_level, a_geom, a_ba, a_dm, a_ngrow, a_support);
}

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::Level* eb_level,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support)
{
    return std::make_unique<EBFArrayBoxFactory>(*eb_level, eb_level->Geom(),
                                                a_ba, a_dm, a_ngrow, a_support);
}

std::unique_ptr<EBFArrayBoxFactory>
makeEBFabFactory (const EB2::IndexSpace* index_space, const Geometry& a_geom,
                  const BoxArray& a_ba,
                  const DistributionMapping& a_dm,
                  const Vector<int>& a_ngrow, EBSupport a_support)
{
    const EB2::Level& eb_level = index_space->getLevel(a_geom);
    return std::make_unique<EBFArrayBoxFactory>(eb_level, a_geom,
                                                a_ba, a_dm, a_ngrow, a_support);
}

}
