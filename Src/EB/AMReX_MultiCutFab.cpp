
#include <AMReX_MultiCutFab.H>
#include <AMReX_MultiFab.H>

#ifdef AMREX_USE_OMP
#include <omp.h>
#endif

namespace amrex {

namespace detail {
    class CutFabFactory final
        : public DefaultFabFactory<CutFab>
    {
    public:
        CutFabFactory (const FabArray<EBCellFlagFab>* a_cellflags)
            : m_cellflags(a_cellflags)
            {}

        AMREX_NODISCARD
        CutFab* create (const Box& box, int ncomps, const FabInfo& info, int box_index) const override
        {
            if ((*m_cellflags)[box_index].getType() == FabType::singlevalued) {
                return new CutFab(box, ncomps, info.alloc, info.shared, info.arena);
            } else {
                return new CutFab();
            }
        }

        AMREX_NODISCARD
        CutFabFactory* clone () const override { return new CutFabFactory(m_cellflags); }

        AMREX_NODISCARD
        Long nBytes (const Box& box, int ncomps, int box_index) const override
        {
            if ((*m_cellflags)[box_index].getType() == FabType::singlevalued) {
                return box.numPts() * ncomps * Long(sizeof(Real));
            } else {
                return Long(0);
            }
        }

    private:
        const FabArray<EBCellFlagFab>* m_cellflags = nullptr;
    };
}


MultiCutFab::MultiCutFab () = default;

MultiCutFab::MultiCutFab (const BoxArray& ba, const DistributionMapping& dm,
                          int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
    : m_data(ba,dm,ncomp,ngrow,MFInfo(),detail::CutFabFactory(&cellflags)),
      m_cellflags(&cellflags)
{
}

MultiCutFab::~MultiCutFab () = default;

void
MultiCutFab::define (const BoxArray& ba, const DistributionMapping& dm,
                     int ncomp, int ngrow, const FabArray<EBCellFlagFab>& cellflags)
{
    m_data.define(ba,dm,ncomp,ngrow,MFInfo(),detail::CutFabFactory(&cellflags)),
    m_cellflags = &cellflags;
}

const CutFab&
MultiCutFab::operator[] (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

CutFab&
MultiCutFab::operator[] (const MFIter& mfi) noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data[mfi];
}

const CutFab&
MultiCutFab::operator[] (int global_box_index) const noexcept
{
    AMREX_ASSERT(ok(global_box_index));
    return m_data[global_box_index];
}

CutFab&
MultiCutFab::operator[] (int global_box_index) noexcept
{
    AMREX_ASSERT(ok(global_box_index));
    return m_data[global_box_index];
}

Array4<Real const>
MultiCutFab::const_array (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.array(mfi);
}

Array4<Real const>
MultiCutFab::array (const MFIter& mfi) const noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.array(mfi);
}

Array4<Real>
MultiCutFab::array (const MFIter& mfi) noexcept
{
    AMREX_ASSERT(ok(mfi));
    return m_data.array(mfi);
}

bool
MultiCutFab::ok (const MFIter& mfi) const noexcept
{
    return (*m_cellflags)[mfi].getType() == FabType::singlevalued;
}

bool
MultiCutFab::ok (int global_box_index) const noexcept
{
    return (*m_cellflags)[global_box_index].getType() == FabType::singlevalued;
}

void
MultiCutFab::setVal (Real val)
{
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(m_data); mfi.isValid(); ++mfi)
    {
        if (ok(mfi)) {
            Array4<Real> const& a = m_data.array(mfi);
            Box const& b = mfi.fabbox();
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(b, m_data.nComp(), i, j, k, n,
            {
                a(i,j,k,n) = val;
            });
        }
    }
}

void
MultiCutFab::ParallelCopy (const MultiCutFab& src, int scomp, int dcomp, int ncomp, int sng, int dng, const Periodicity& period)
{
    m_data.ParallelCopy(src.m_data, scomp, dcomp, ncomp, sng, dng, period);
}

MultiFab
MultiCutFab::ToMultiFab (Real regular_value, Real covered_value) const
{
    const int ncomp = nComp();
    MultiFab mf(boxArray(), DistributionMap(), ncomp, nGrow());
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        auto t = (*m_cellflags)[mfi].getType();
        Box const& b = mfi.fabbox();
        Array4<Real> const& d = mf.array(mfi);
        if (t == FabType::singlevalued) {
            Array4<Real const> const& s = m_data.const_array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(b, ncomp, i, j, k, n,
            {
                d(i,j,k,n) = s(i,j,k,n);
            });
        } else {
            Real val = (t == FabType::regular) ? regular_value : covered_value;
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D(b, ncomp, i, j, k, n,
            {
                d(i,j,k,n) = val;
            });
        }
    }
    return mf;
}

}
