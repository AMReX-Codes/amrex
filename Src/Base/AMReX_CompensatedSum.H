#ifndef AMREX_COMPENSATEDSUM_H_
#define AMREX_COMPENSATEDSUM_H_
#include <AMReX_Config.H>

#include <AMReX_BLassert.H>
#include <AMReX_GpuQualifiers.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>

#include <type_traits>
#include <utility>

namespace amrex
{

template <typename T, typename AccumT = T>
struct CompensatedSum
{
    AccumT sum = AccumT(0);
    AccumT comp = AccumT(0);

    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (U value) noexcept
    {
        AccumT const y = static_cast<AccumT>(value) - comp;
        AccumT const t = sum + y;
        comp = (t - sum) - y;
        sum = t;
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (CompensatedSum const& rhs) noexcept
    {
        add(rhs.value());
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    AccumT value () const noexcept { return sum; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    AccumT compensation () const noexcept { return comp; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void reset () noexcept
    {
        sum = AccumT(0);
        comp = AccumT(0);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (T value) noexcept
    {
        sum = static_cast<AccumT>(value);
        comp = AccumT(0);
    }
};

template <typename T, typename U>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void compensatedAdd (T& sum, T& comp, U value) noexcept
{
    T const y = static_cast<T>(value) - comp;
    T const t = sum + y;
    comp = (t - sum) - y;
    sum = t;
}

template <typename T>
struct CompensatedRef
{
    T* sum = nullptr;
    T* comp = nullptr;

    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (U value) const noexcept
    {
        compensatedAdd(*sum, *comp, value);
    }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T value () const noexcept { return *sum; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    T compensation () const noexcept { return *comp; }

    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void set (T value) const noexcept
    {
        *sum = value;
        *comp = T(0);
    }
};

template <typename SumView, typename CompView>
struct CompensatedView
{
    SumView sum;
    CompView comp;

    template <typename... I>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    auto operator() (I... idx) const noexcept
    {
        using T = std::remove_cv_t<std::remove_reference_t<decltype(sum(idx...))>>;
        return CompensatedRef<T>{&(sum(idx...)), &(comp(idx...))};
    }

    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (int i, int j, int k, U value) const noexcept
    {
        (*this)(i,j,k).add(value);
    }

    template <typename U>
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    void add (int i, int j, int k, int n, U value) const noexcept
    {
        (*this)(i,j,k,n).add(value);
    }
};

template <typename SumView, typename CompView>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto makeCompensatedView (SumView const& sum, CompView const& comp) noexcept
{
    return CompensatedView<SumView,CompView>{sum, comp};
}

class CompensatedMultiFab
{
public:
    CompensatedMultiFab () noexcept = default;

    CompensatedMultiFab (CompensatedMultiFab const&) = delete;
    CompensatedMultiFab& operator= (CompensatedMultiFab const&) = delete;

    CompensatedMultiFab (CompensatedMultiFab&& rhs)
        : m_owned_value(std::move(rhs.m_owned_value)),
          m_owned_comp(std::move(rhs.m_owned_comp))
    {
        m_value = (rhs.m_value == &rhs.m_owned_value) ? &m_owned_value : rhs.m_value;
        m_comp = (rhs.m_comp == &rhs.m_owned_comp) ? &m_owned_comp : rhs.m_comp;
        rhs.m_value = nullptr;
        rhs.m_comp = nullptr;
    }

    CompensatedMultiFab& operator= (CompensatedMultiFab&& rhs)
    {
        if (this != &rhs) {
            bool const value_is_owned = (rhs.m_value == &rhs.m_owned_value);
            bool const comp_is_owned = (rhs.m_comp == &rhs.m_owned_comp);
            m_owned_value = std::move(rhs.m_owned_value);
            m_owned_comp = std::move(rhs.m_owned_comp);
            m_value = value_is_owned ? &m_owned_value : rhs.m_value;
            m_comp = comp_is_owned ? &m_owned_comp : rhs.m_comp;
            rhs.m_value = nullptr;
            rhs.m_comp = nullptr;
        }
        return *this;
    }

    explicit CompensatedMultiFab (MultiFab& value)
    {
        define(value);
    }

    CompensatedMultiFab (MultiFab& value, MultiFab& compensation) noexcept
    {
        define(value, compensation);
    }

    CompensatedMultiFab (BoxArray const& bxs, DistributionMapping const& dm,
                         int ncomp, int ngrow,
                         MFInfo const& info = MFInfo(),
                         FabFactory<FArrayBox> const& factory = FArrayBoxFactory())
    {
        define(bxs, dm, ncomp, ngrow, info, factory);
    }

    CompensatedMultiFab (BoxArray const& bxs, DistributionMapping const& dm,
                         int ncomp, IntVect const& ngrow,
                         MFInfo const& info = MFInfo(),
                         FabFactory<FArrayBox> const& factory = FArrayBoxFactory())
    {
        define(bxs, dm, ncomp, ngrow, info, factory);
    }

    void define (MultiFab& value)
    {
        m_owned_value = MultiFab();
        m_owned_comp.define(value.boxArray(), value.DistributionMap(), value.nComp(),
                            value.nGrowVect(), MFInfo().SetArena(value.arena()),
                            value.Factory());
        m_value = &value;
        m_comp = &m_owned_comp;
        resetCompensation();
    }

    void define (MultiFab& value, MultiFab& compensation) noexcept
    {
        BL_ASSERT(value.boxArray() == compensation.boxArray());
        BL_ASSERT(value.DistributionMap() == compensation.DistributionMap());
        BL_ASSERT(value.nComp() == compensation.nComp());
        BL_ASSERT(value.nGrowVect() == compensation.nGrowVect());
        m_owned_value = MultiFab();
        m_owned_comp = MultiFab();
        m_value = &value;
        m_comp = &compensation;
    }

    void define (BoxArray const& bxs, DistributionMapping const& dm,
                 int ncomp, int ngrow,
                 MFInfo const& info = MFInfo(),
                 FabFactory<FArrayBox> const& factory = FArrayBoxFactory())
    {
        m_owned_value.define(bxs, dm, ncomp, ngrow, info, factory);
        m_owned_comp.define(bxs, dm, ncomp, ngrow, info, factory);
        m_value = &m_owned_value;
        m_comp = &m_owned_comp;
        setVal(Real(0.0));
    }

    void define (BoxArray const& bxs, DistributionMapping const& dm,
                 int ncomp, IntVect const& ngrow,
                 MFInfo const& info = MFInfo(),
                 FabFactory<FArrayBox> const& factory = FArrayBoxFactory())
    {
        m_owned_value.define(bxs, dm, ncomp, ngrow, info, factory);
        m_owned_comp.define(bxs, dm, ncomp, ngrow, info, factory);
        m_value = &m_owned_value;
        m_comp = &m_owned_comp;
        setVal(Real(0.0));
    }

    void resetCompensation ()
    {
        BL_ASSERT(m_comp != nullptr);
        m_comp->setVal(Real(0.0));
    }

    void setVal (Real value)
    {
        BL_ASSERT(m_value != nullptr);
        m_value->setVal(value);
        resetCompensation();
    }

    [[nodiscard]] MultiFab& value () noexcept
    {
        BL_ASSERT(m_value != nullptr);
        return *m_value;
    }

    [[nodiscard]] MultiFab const& value () const noexcept
    {
        BL_ASSERT(m_value != nullptr);
        return *m_value;
    }

    [[nodiscard]] MultiFab& compensation () noexcept
    {
        BL_ASSERT(m_comp != nullptr);
        return *m_comp;
    }

    [[nodiscard]] MultiFab const& compensation () const noexcept
    {
        BL_ASSERT(m_comp != nullptr);
        return *m_comp;
    }

    [[nodiscard]] auto array (MFIter const& mfi) noexcept
    {
        return makeCompensatedView(m_value->array(mfi), m_comp->array(mfi));
    }

    [[nodiscard]] auto array (MFIter const& mfi, int start_comp) noexcept
    {
        return makeCompensatedView(m_value->array(mfi, start_comp),
                                   m_comp->array(mfi, start_comp));
    }

    [[nodiscard]] auto array (MFIter const& mfi, int start_comp, int num_comps) noexcept
    {
        return makeCompensatedView((*m_value)[mfi].array(start_comp, num_comps),
                                   (*m_comp)[mfi].array(start_comp, num_comps));
    }

private:
    MultiFab m_owned_value;
    MultiFab m_owned_comp;
    MultiFab* m_value = nullptr;
    MultiFab* m_comp = nullptr;
};

}

#endif
