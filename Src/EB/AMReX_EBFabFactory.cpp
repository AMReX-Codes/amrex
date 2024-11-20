
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
{
    auto const& ebflags = getMultiEBCellFlagFab();
#ifdef AMREX_USE_GPU
    m_eb_data.resize(EBData::real_data_size*ebflags.local_size());
    Gpu::PinnedVector<Array4<Real const>> eb_data_hv;
#else
    auto& eb_data_hv = m_eb_data;
#endif

    eb_data_hv.reserve(EBData::real_data_size*ebflags.local_size());

    for (MFIter mfi(ebflags,MFItInfo{}.DisableDeviceSync()); mfi.isValid(); ++mfi) {
        Array4<Real const> a{};

        bool cutfab_is_ok = ebflags[mfi].getType() == FabType::singlevalued;

        a = ( m_ebdc->m_levelset )
            ? m_ebdc->m_levelset->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        a = ( m_ebdc->m_volfrac )
            ? m_ebdc->m_volfrac->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        a = ( m_ebdc->m_centroid && cutfab_is_ok )
            ? m_ebdc->m_centroid->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        a = ( m_ebdc->m_bndrycent && cutfab_is_ok )
            ? m_ebdc->m_bndrycent->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        a = ( m_ebdc->m_bndrynorm && cutfab_is_ok )
            ? m_ebdc->m_bndrynorm->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        a = ( m_ebdc->m_bndryarea && cutfab_is_ok )
            ? m_ebdc->m_bndryarea->const_array(mfi) : Array4<Real const>{};
        eb_data_hv.push_back(a);

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            a = ( m_ebdc->m_areafrac[idim] && cutfab_is_ok )
                ? m_ebdc->m_areafrac[idim]->const_array(mfi) : Array4<Real const>{};
            eb_data_hv.push_back(a);
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            a = ( m_ebdc->m_facecent[idim] && cutfab_is_ok )
                ? m_ebdc->m_facecent[idim]->const_array(mfi) : Array4<Real const>{};
            eb_data_hv.push_back(a);
        }

        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            a = ( m_ebdc->m_edgecent[idim] && cutfab_is_ok )
                ? m_ebdc->m_edgecent[idim]->const_array(mfi) : Array4<Real const>{};
            eb_data_hv.push_back(a);
        }
    }

#ifdef AMREX_USE_GPU
    Gpu::copyAsync(Gpu::hostToDevice, eb_data_hv.begin(), eb_data_hv.end(), m_eb_data.begin());
    Gpu::streamSynchronize();
#endif
}

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
        return new EBFArrayBox(ebcellflag, box, ncomps, info.arena, this, box_index);
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
        auto const& ebrhs = static_cast<EBFArrayBox const&>(rhs);
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
        auto* p = static_cast<EBFArrayBox*>(fab);
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

bool
EBFArrayBoxFactory::hasEBInfo () const noexcept
{
    return m_parent->hasEBInfo();
}

EBData
EBFArrayBoxFactory::getEBData (MFIter const& mfi) const noexcept
{
    int const li = mfi.LocalIndex();
    auto const& ebflags_ma = this->getMultiEBCellFlagFab().const_arrays();
#ifdef AMREX_USE_GPU
    auto const* pebflag = ebflags_ma.dp + li;
#else
    auto const* pebflag = ebflags_ma.hp + li;
#endif
    return EBData{pebflag, m_eb_data.data()+EBData::real_data_size*li};
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
