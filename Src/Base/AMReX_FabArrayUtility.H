#ifndef AMREX_FABARRAY_UTILITY_H_
#define AMREX_FABARRAY_UTILITY_H_
#include <AMReX_Config.H>

#include <AMReX_FabArray.H>
#include <AMReX_LayoutData.H>
#include <AMReX_Print.H>
#include <AMReX_ParReduce.H>
#include <limits>

namespace amrex {

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, int nghost, F&& f) {
    return ReduceSum(fa, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum_host (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        auto const& arr = fa.const_array(mfi);
        sm += f(bx, arr);
    }

    return sm;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class OP, class FAB, class F>
std::enable_if_t<IsBaseFab<FAB>::value,
                 std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                    std::is_same<OP,ReduceOpLogicalOr>::value,
                                    int, typename FAB::value_type> >
ReduceMF (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    using T = std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                 std::is_same<OP,ReduceOpLogicalOr>::value,
                                 int, typename FAB::value_type>;
    auto typ = fa.ixType();
    auto const& ma = fa.const_arrays();
    return ParReduce(TypeList<OP>{}, TypeList<T>{}, fa, nghost,
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
           {
               return { static_cast<T>(f(amrex::makeSingleCellBox(i,j,k,typ), ma[box_no])) };
           });
}

template <class OP, class FAB1, class FAB2, class F>
std::enable_if_t<IsBaseFab<FAB1>::value && IsBaseFab<FAB2>::value,
                 std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                    std::is_same<OP,ReduceOpLogicalOr>::value,
                                    int, typename FAB1::value_type> >
ReduceMF (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, IntVect const& nghost, F const& f)
{
    using T = std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                 std::is_same<OP,ReduceOpLogicalOr>::value,
                                 int, typename FAB1::value_type>;
    auto typ = fa1.ixType();
    auto const& ma1 = fa1.const_arrays();
    auto const& ma2 = fa2.const_arrays();
    return ParReduce(TypeList<OP>{}, TypeList<T>{}, fa1, nghost,
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
           {
               return { static_cast<T>(f(amrex::makeSingleCellBox(i,j,k,typ),
                                         ma1[box_no], ma2[box_no])) };
           });
}

template <class OP, class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<IsBaseFab<FAB1>::value && IsBaseFab<FAB2>::value && IsBaseFab<FAB3>::value,
                 std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                    std::is_same<OP,ReduceOpLogicalOr>::value,
                                    int, typename FAB1::value_type> >
ReduceMF (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
          FabArray<FAB3> const& fa3, IntVect const& nghost, F const& f)
{
    using T = std::conditional_t<std::is_same<OP,ReduceOpLogicalAnd>::value ||
                                 std::is_same<OP,ReduceOpLogicalOr>::value,
                                 int, typename FAB1::value_type>;
    auto typ = fa1.ixType();
    auto const& ma1 = fa1.const_arrays();
    auto const& ma2 = fa2.const_arrays();
    auto const& ma3 = fa3.const_arrays();
    return ParReduce(TypeList<OP>{}, TypeList<T>{}, fa1, nghost,
           [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
           {
               return { static_cast<T>(f(amrex::makeSingleCellBox(i,j,k,typ),
                                         ma1[box_no], ma2[box_no], ma3[box_no])) };
           });
}

template <class FAB, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceSum_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return ReduceSum_host(fa,nghost,std::forward<F>(f));
}

template <class FAB, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceSum_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa,nghost,f);
    amrex::Abort("ReduceSum: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpSum>(fa, nghost, std::forward<F>(f));
    } else {
        return fudetail::ReduceSum_host_wrapper(fa, nghost, std::forward<F>(f));
    }
}
#else
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceSum (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceSum_host(fa, nghost, std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           int nghost, F&& f) {
    return ReduceSum(fa1, fa2, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        sm += f(bx, arr1, arr2);
    }

    return sm;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceSum_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    return ReduceSum_host(fa1,fa2,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceSum_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,nghost,f);
    amrex::Abort("ReduceSum: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpSum>(fa1,fa2,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceSum_host_wrapper(fa1,fa2,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    return fudetail::ReduceSum_host(fa1,fa2,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, FabArray<FAB3> const& fa3,
           int nghost, F&& f)
{
  return ReduceSum(fa1, fa2, fa3, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                FabArray<FAB3> const& fa3, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type sm = 0;

#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        const auto& arr3 = fa3.const_array(mfi);
        sm += f(bx, arr1, arr2, arr3);
    }

    return sm;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceSum_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceSum_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceSum_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,fa3,nghost,f);
    amrex::Abort("ReduceSum: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpSum>(fa1,fa2,fa3,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceSum_host_wrapper(fa1,fa2,fa3,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceSum (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceSum_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}
#endif

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin (FabArray<FAB> const& fa, int nghost, F&& f)
{
    return ReduceMin(fa, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin_host (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB::value_type;
    value_type r = std::numeric_limits<value_type>::max();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:r)
#endif
    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr = fa.const_array(mfi);
        r = std::min(r, f(bx, arr));
    }
    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceMin_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return ReduceMin_host(fa,nghost,std::forward<F>(f));
}

template <class FAB, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceMin_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa,nghost,f);
    amrex::Abort("ReduceMin: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMin>(fa, nghost, std::forward<F>(f));
    } else {
        return fudetail::ReduceMin_host_wrapper(fa, nghost, std::forward<F>(f));
    }
}
#else
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMin (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMin_host(fa, nghost, std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, int nghost, F&& f)
{
    return ReduceMin(fa1, fa2, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::max();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        r = std::min(r, f(bx, arr1, arr2));
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMin_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMin_host(fa1,fa2,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMin_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,nghost,f);
    amrex::Abort("ReduceMin: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMin>(fa1,fa2,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceMin_host_wrapper(fa1,fa2,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMin_host(fa1,fa2,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, FabArray<FAB3> const& fa3,
           int nghost, F&& f)
{
    return ReduceMin(fa1, fa2, fa3, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                FabArray<FAB3> const& fa3, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::max();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(min:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        const auto& arr3 = fa3.const_array(mfi);
        r = std::min(r, f(bx, arr1, arr2, arr3));
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMin_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMin_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMin_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,fa3,nghost,f);
    amrex::Abort("ReduceMin: Launch Region is off. Device lambda lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMin>(fa1,fa2,fa3,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceMin_host_wrapper(fa1,fa2,fa3,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMin (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMin_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}
#endif

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax (FabArray<FAB> const& fa, int nghost, F&& f)
{
    return ReduceMax(fa, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax_host (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB::value_type;
    value_type r = std::numeric_limits<value_type>::lowest();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:r)
#endif
    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr = fa.const_array(mfi);
        r = std::max(r, f(bx, arr));
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceMax_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return ReduceMax_host(fa,nghost,std::forward<F>(f));
}

template <class FAB, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB::value_type>
ReduceMax_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa,nghost,f);
    amrex::Abort("ReduceMax: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMax>(fa,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceMax_host_wrapper(fa,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
typename FAB::value_type
ReduceMax (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMax_host(fa,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, int nghost, F&& f)
{
    return ReduceMax(fa1, fa2, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::lowest();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        r = std::max(r, f(bx, arr1, arr2));
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMax_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    return ReduceMax_host(fa1,fa2,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMax_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,nghost,f);
    amrex::Abort("ReduceMax: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMax>(fa1,fa2,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceMax_host_wrapper(fa1,fa2,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMax_host(fa1,fa2,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2, FabArray<FAB3> const& fa3,
           int nghost, F&& f)
{
    return ReduceMax(fa1, fa2, fa3, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                FabArray<FAB3> const& fa3, IntVect const& nghost, F const& f)
{
    using value_type = typename FAB1::value_type;
    value_type r = std::numeric_limits<value_type>::lowest();

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(max:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        const auto& arr3 = fa3.const_array(mfi);
        r = std::max(r, f(bx, arr1, arr2, arr3));
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMax_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMax_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class FAB3, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, typename FAB1::value_type>
ReduceMax_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                        FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,fa3,nghost,f);
    amrex::Abort("ReduceMax: Launch Region is off. Device lambda lambda cannot be called by host.");
    return 0;
}
}

template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpMax>(fa1,fa2,fa3,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceMax_host_wrapper(fa1,fa2,fa3,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class FAB3, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
typename FAB1::value_type
ReduceMax (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
           FabArray<FAB3> const& fa3, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceMax_host(fa1,fa2,fa3,nghost,std::forward<F>(f));
}
#endif

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd (FabArray<FAB> const& fa, int nghost, F&& f)
{
    return ReduceLogicalAnd(fa, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd_host (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    int r = true;

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(&&:r)
#endif
    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr = fa.const_array(mfi);
        r = r && f(bx, arr);
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalAnd_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return ReduceLogicalAnd_host(fa,nghost,std::forward<F>(f));
}

template <class FAB, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalAnd_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa,nghost,f);
    amrex::Abort("ReduceLogicalAnd: Launch Region is off. Device lambda cannot be called by host.");
    return false;
}
}

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpLogicalAnd>(fa,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceLogicalAnd_host_wrapper(fa,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalAnd (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceLogicalAnd_host(fa,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  int nghost, F&& f)
{
    return ReduceLogicalAnd(fa1, fa2, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                       IntVect const& nghost, F const& f)
{
    int r = true;

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(&&:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        r = r && f(bx, arr1, arr2);
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalAnd_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  IntVect const& nghost, F&& f)
{
    return ReduceLogicalAnd_host(fa1,fa2,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalAnd_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                               IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,nghost,f);
    amrex::Abort("ReduceLogicalAnd: Luanch Region is off. Device lambda cannot be called by host.");
    return false;
}
}

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpLogicalAnd>(fa1,fa2,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceLogicalAnd_host_wrapper(fa1,fa2,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalAnd (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                  IntVect const& nghost, F&& f)
{
    return fudetail::ReduceLogicalAnd_host(fa1,fa2,nghost,std::forward<F>(f));
}
#endif

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr (FabArray<FAB> const& fa, int nghost, F&& f)
{
    return ReduceLogicalOr(fa, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr_host (FabArray<FAB> const& fa, IntVect const& nghost, F const& f)
{
    int r = false;

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(||:r)
#endif
    for (MFIter mfi(fa,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr = fa.const_array(mfi);
        r = r || f(bx, arr);
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalOr_host_wrapper (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return ReduceLogicalOr_host(fa,nghost,std::forward<F>(f));
}

template <class FAB, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalOr_host (FabArray<FAB> const& fa, IntVect const& nghost, F&& /*f*/)
{
    amrex::ignore_unused(fa,nghost);
    amrex::Abort("ReduceLogicalOr: Launch Region is off. Device lambda cannot be called by host.");
    return 0;
}
}

template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpLogicalOr>(fa,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceLogicalOr_host_wrapper(fa,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
bool
ReduceLogicalOr (FabArray<FAB> const& fa, IntVect const& nghost, F&& f)
{
    return fudetail::ReduceLogicalOr_host(fa,nghost,std::forward<F>(f));
}
#endif

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 int nghost, F&& f)
{
    return ReduceLogicalOr(fa1, fa2, IntVect(nghost), std::forward<F>(f));
}

namespace fudetail {
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr_host (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                      IntVect const& nghost, F const& f)
{
    int r = false;

#ifdef AMREX_USE_OMP
#pragma omp parallel reduction(||:r)
#endif
    for (MFIter mfi(fa1,true); mfi.isValid(); ++mfi)
    {
        const Box& bx = mfi.growntilebox(nghost);
        const auto& arr1 = fa1.const_array(mfi);
        const auto& arr2 = fa2.const_array(mfi);
        r = r || f(bx, arr1, arr2);
    }

    return r;
}
}

#ifdef AMREX_USE_GPU
namespace fudetail {
template <class FAB1, class FAB2, class F>
std::enable_if_t<!amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalOr_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 IntVect const& nghost, F&& f)
{
    return fudetail::ReduceLogicalOr_host(fa1,fa2,nghost,std::forward<F>(f));
}

template <class FAB1, class FAB2, class F>
std::enable_if_t<amrex::DefinitelyNotHostRunnable<F>::value, bool>
ReduceLogicalOr_host_wrapper (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                              IntVect const& nghost, F&& f)
{
    amrex::ignore_unused(fa1,fa2,nghost,f);
    amrex::Abort("ReeuceLogicalOr: Launch Region is off. Device lambda cannot be called by host.");
    return false;
}
}

template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 IntVect const& nghost, F&& f)
{
    if (Gpu::inLaunchRegion()) {
        return fudetail::ReduceMF<ReduceOpLogicalOr>(fa1,fa2,nghost,std::forward<F>(f));
    } else {
        return fudetail::ReduceLogicalOr_host_wrapper(fa1,fa2,nghost,std::forward<F>(f));
    }
}
#else
template <class FAB1, class FAB2, class F,
          class bar = std::enable_if_t<IsBaseFab<FAB1>::value> >
bool
ReduceLogicalOr (FabArray<FAB1> const& fa1, FabArray<FAB2> const& fa2,
                 IntVect const& nghost, F&& f)
{
    return fudetail::ReduceLogicalOr_host(fa1,fa2,nghost,std::forward<F>(f));
}
#endif

template <class FAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
printCell (FabArray<FAB> const& mf, const IntVect& cell, int comp = -1,
           const IntVect& ng = IntVect::TheZeroVector())
{
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box& bx = amrex::grow(mfi.validbox(), ng);
        if (bx.contains(cell)) {
            int n = (comp >= 0) ? 1 : mf.nComp();
            auto const& fab = mf.const_array(mfi);
            Gpu::PinnedVector<typename FAB::value_type> pv(n);
            auto* dp = pv.data();
            auto f = [=] AMREX_GPU_HOST_DEVICE ()
                {
                    if (comp >= 0) {
                        *dp = fab(cell, comp);
                    } else {
                        for (int i = 0; i < n; ++i) {
                            dp[i] = fab(cell,i);
                        }
                    }
                };

#ifdef AMREX_USE_GPU
            if (mf.arena()->isManaged() || mf.arena()->isDevice()) {
                amrex::single_task(f);
                Gpu::streamSynchronize();
            } else
#endif
            {
                f();
            }

            if (comp >= 0) {
                amrex::AllPrint().SetPrecision(17) << " At cell " << cell << " in Box " << bx
                                                   << ": " << *dp << '\n';
            } else {
                std::ostringstream ss;
                ss.precision(17);
                for (int i = 0; i < n-1; ++i)
                {
                    ss << dp[i] << ", ";
                }
                ss << dp[n-1];
                amrex::AllPrint() << " At cell " << cell << " in Box " << bx
                                  << ": " << ss.str() << '\n';
            }
        }
    }
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Subtract (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Subtract(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Subtract (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstfa = dst.arrays();
        auto const& srcfa = src.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstfa[box_no](i,j,k,n+dstcomp) -= srcfa[box_no](i,j,k,n+srccomp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const srcFab = src.array(mfi);
                auto       dstFab = dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dstFab(i,j,k,n+dstcomp) -= srcFab(i,j,k,n+srccomp);
                });
            }
        }
    }
}


template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Multiply (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Multiply(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Multiply (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstfa = dst.arrays();
        auto const& srcfa = src.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstfa[box_no](i,j,k,n+dstcomp) *= srcfa[box_no](i,j,k,n+srccomp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const srcFab = src.array(mfi);
                auto       dstFab = dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dstFab(i,j,k,n+dstcomp) *= srcFab(i,j,k,n+srccomp);
                });
            }
        }
    }
}


template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Divide (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, int nghost)
{
    Divide(dst,src,srccomp,dstcomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Divide (FabArray<FAB>& dst, FabArray<FAB> const& src, int srccomp, int dstcomp, int numcomp, const IntVect& nghost)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && dst.isFusingCandidate()) {
        auto const& dstfa = dst.arrays();
        auto const& srcfa = src.const_arrays();
        ParallelFor(dst, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            dstfa[box_no](i,j,k,n+dstcomp) /= srcfa[box_no](i,j,k,n+srccomp);
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(dst,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const srcFab = src.array(mfi);
                auto       dstFab = dst.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    dstFab(i,j,k,n+dstcomp) /= srcFab(i,j,k,n+srccomp);
                });
            }
        }
    }
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Abs (FabArray<FAB>& fa, int icomp, int numcomp, int nghost)
{
    Abs(fa,icomp,numcomp,IntVect(nghost));
}

template <class FAB,
          class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
Abs (FabArray<FAB>& fa, int icomp, int numcomp, const IntVect& nghost)
{
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && fa.isFusingCandidate()) {
        auto const& fabarr = fa.arrays();
        ParallelFor(fa, nghost, numcomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            fabarr[box_no](i,j,k,n+icomp) = std::abs(fabarr[box_no](i,j,k,n+icomp));
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(fa,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.growntilebox(nghost);
            if (bx.ok())
            {
                auto const& fab = fa.array(mfi);
                AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, numcomp, i, j, k, n,
                {
                    fab(i,j,k,n+icomp) = std::abs(fab(i,j,k,n+icomp));
                });
            }
        }
    }
}

template <class FAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
prefetchToHost (FabArray<FAB> const& fa, const bool synchronous = true)
{
#ifdef AMREX_USE_GPU
    if (fa.arena()->isManaged()) {
        for (MFIter mfi(fa, MFItInfo().SetDeviceSync(synchronous)); mfi.isValid(); ++mfi) {
            fa.prefetchToHost(mfi);
        }
    }
#else
    amrex::ignore_unused(fa,synchronous);
#endif
}

template <class FAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
prefetchToDevice (FabArray<FAB> const& fa, const bool synchronous = true)
{
#ifdef AMREX_USE_GPU
    if (fa.arena()->isManaged()) {
        for (MFIter mfi(fa, MFItInfo().SetDeviceSync(synchronous)); mfi.isValid(); ++mfi) {
            fa.prefetchToDevice(mfi);
        }
    }
#else
    amrex::ignore_unused(fa,synchronous);
#endif
}


template <class FAB, class IFAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value
                                                               && IsBaseFab<IFAB>::value> >
void
OverrideSync (FabArray<FAB> & fa, FabArray<IFAB> const& msk, const Periodicity& period)
{
    BL_PROFILE("OverrideSync()");

    OverrideSync_nowait(fa, msk, period);
    OverrideSync_finish(fa);
}


template <class FAB, class IFAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value
                                                               && IsBaseFab<IFAB>::value> >
void
OverrideSync_nowait (FabArray<FAB> & fa, FabArray<IFAB> const& msk, const Periodicity& period)
{
    BL_PROFILE("OverrideSync_nowait()");
    AMREX_ASSERT_WITH_MESSAGE(!fa.os_temp, "OverrideSync_nowait() called when already in progress.");

    if (fa.ixType().cellCentered()) { return; }

    const int ncomp = fa.nComp();

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion() && fa.isFusingCandidate()) {
        auto const&  fabarr = fa.arrays();
        auto const& ifabarr = msk.const_arrays();
        ParallelFor(fa, IntVect(0), ncomp,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k, int n) noexcept
        {
            if (!ifabarr[box_no](i,j,k)) { fabarr[box_no](i,j,k,n) = 0; }
        });
        if (!Gpu::inNoSyncRegion()) {
            Gpu::streamSynchronize();
        }
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        for (MFIter mfi(fa,TilingIfNotGPU()); mfi.isValid(); ++mfi)
        {
            const Box& bx = mfi.tilebox();
            auto fab = fa.array(mfi);
            auto const ifab = msk.array(mfi);
            AMREX_HOST_DEVICE_PARALLEL_FOR_4D( bx, ncomp, i, j, k, n,
            {
                if (!ifab(i,j,k)) { fab(i,j,k,n) = 0; }
            });
        }
    }

    fa.os_temp = std::make_unique< FabArray<FAB> > ( fa.boxArray(), fa.DistributionMap(),
                                                     ncomp, 0, MFInfo(), fa.Factory() );
    fa.os_temp->setVal(0);
    fa.os_temp->ParallelCopy_nowait(fa, period, FabArrayBase::ADD);
}

template <class FAB, class bar = std::enable_if_t<IsBaseFab<FAB>::value> >
void
OverrideSync_finish (FabArray<FAB> & fa)
{
    BL_PROFILE("OverrideSync_finish()");

    if (fa.ixType().cellCentered()) { return; }

    fa.os_temp->ParallelCopy_finish();
    amrex::Copy(fa, *(fa.os_temp), 0, 0, fa.nComp(), 0);

    fa.os_temp.reset();
}

template <class FAB, class foo = std::enable_if_t<IsBaseFab<FAB>::value> >
void
dtoh_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src,
             int scomp, int dcomp, int ncomp)
{
    AMREX_ASSERT(isMFIterSafe(dst, src));
    AMREX_ASSERT(dst.nGrowVect() == src.nGrowVect());
#ifdef AMREX_USE_GPU
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        void* pdst = dst[mfi].dataPtr(dcomp);
        void const* psrc = src[mfi].dataPtr(scomp);
        Gpu::dtoh_memcpy_async(pdst, psrc, dst[mfi].nBytes(mfi.fabbox(), ncomp));
    }
#else
    Copy(dst, src, scomp, dcomp, ncomp, dst.nGrowVect());
#endif
}

template <class FAB, class foo = std::enable_if_t<IsBaseFab<FAB>::value> >
void
dtoh_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src)
{
    dtoh_memcpy(dst, src, 0, 0, dst.nComp());
}

template <class FAB, class foo = std::enable_if_t<IsBaseFab<FAB>::value> >
void
htod_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src,
             int scomp, int dcomp, int ncomp)
{
    AMREX_ASSERT(isMFIterSafe(dst, src));
    AMREX_ASSERT(dst.nGrowVect() == src.nGrowVect());
#ifdef AMREX_USE_GPU
    for (MFIter mfi(dst); mfi.isValid(); ++mfi) {
        void* pdst = dst[mfi].dataPtr(dcomp);
        void const* psrc = src[mfi].dataPtr(scomp);
        Gpu::htod_memcpy_async(pdst, psrc, dst[mfi].nBytes(mfi.fabbox(), ncomp));
    }
#else
    Copy(dst, src, scomp, dcomp, ncomp, dst.nGrowVect());
#endif
}

template <class FAB, class foo = std::enable_if_t<IsBaseFab<FAB>::value> >
void
htod_memcpy (FabArray<FAB>& dst, FabArray<FAB> const& src)
{
    htod_memcpy(dst, src, 0, 0, dst.nComp());
}

template <class FAB, class foo = std::enable_if_t<IsBaseFab<FAB>::value> >
IntVect
indexFromValue (FabArray<FAB> const& mf, int comp, IntVect const& nghost,
                typename FAB::value_type value)
{
    IntVect loc;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        amrex::Gpu::Buffer<int> aa({0,AMREX_D_DECL(0,0,0)});
        int* p = aa.data();
        // This is a device ptr to 1+AMREX_SPACEDIM int zeros.
        // The first is used as an atomic bool and the others for intvect.
        if (mf.isFusingCandidate()) {
            auto const& ma = mf.const_arrays();
            ParallelFor(mf, nghost, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
            {
                int* flag = p;
                if (*flag == 0) {
                    if (ma[box_no](i,j,k,comp) == value) {
                        if (Gpu::Atomic::Exch(flag,1) == 0) {
                            AMREX_D_TERM(p[1] = i;,
                                         p[2] = j;,
                                         p[3] = k;);
                        }
                    }
                }
            });
        } else {
            for (MFIter mfi(mf,MFItInfo().SetDeviceSync(false)); mfi.isValid(); ++mfi) {
                const Box& bx = amrex::grow(mfi.validbox(), nghost);
                auto const& arr = mf.const_array(mfi);
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    int* flag = p;
                    if (*flag == 0) {
                        if (arr(i,j,k,comp) == value) {
                            if (Gpu::Atomic::Exch(flag,1) == 0) {
                                AMREX_D_TERM(p[1] = i;,
                                             p[2] = j;,
                                             p[3] = k;);
                            }
                        }
                    }
                });
            }
        }
        int const* tmp = aa.copyToHost();
        AMREX_D_TERM(loc[0] = tmp[1];,
                     loc[1] = tmp[2];,
                     loc[2] = tmp[3];);
    }
    else
#endif
    {
        bool f = false;
#ifdef AMREX_USE_OMP
#pragma omp parallel
#endif
        {
            IntVect priv_loc = IntVect::TheMinVector();
            for (MFIter mfi(mf,true); mfi.isValid(); ++mfi)
            {
                const Box& bx = mfi.growntilebox(nghost);
                auto const& fab = mf.const_array(mfi);
                AMREX_LOOP_3D(bx, i, j, k,
                {
                    if (fab(i,j,k,comp) == value) {
                        priv_loc = IntVect(AMREX_D_DECL(i,j,k));
                    }
                });
            }

            if (priv_loc.allGT(IntVect::TheMinVector())) {
                bool old;
// we should be able to test on _OPENMP < 201107 for capture (version 3.1)
// but we must work around a bug in gcc < 4.9
// And, with NVHPC 21.9 to <23.1, we saw an ICE with the atomic capture (NV bug: #3390723)
#if defined(AMREX_USE_OMP) && defined(_OPENMP) && (_OPENMP < 201307 || (defined(__NVCOMPILER) && __NVCOMPILER_MAJOR__ < 23)) // OpenMP 4.0
#pragma omp critical (amrex_indexfromvalue)
#elif defined(AMREX_USE_OMP)
#pragma omp atomic capture
#endif
                {
                    old = f;
                    f = true;
                }

                if (old == false) { loc = priv_loc; }
            }
        }
    }

    return loc;
}

/**
 * \brief Compute dot products of two FabArrays
 *
 * \param x      first FabArray
 * \param xcomp  starting component of x
 * \param y      second FabArray
 * \param ycomp  starting component of y
 * \param ncomp  number of components
 * \param nghost number of ghost cells
 * \param local  If true, MPI communication is skipped.
 */
template <typename FAB, std::enable_if_t<IsBaseFab<FAB>::value,int> FOO = 0>
typename FAB::value_type
Dot (FabArray<FAB> const& x, int xcomp, FabArray<FAB> const& y, int ycomp, int ncomp,
     IntVect const& nghost, bool local = false)
{
    BL_ASSERT(x.boxArray() == y.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap());
    BL_ASSERT(x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost));

    BL_PROFILE("amrex::Dot()");

    using T = typename FAB::value_type;
    auto sm = T(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<T>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
        {
            auto t = T(0.0);
            auto const& xfab = xma[box_no];
            auto const& yfab = yma[box_no];
            for (int n = 0; n < ncomp; ++n) {
                t += xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& xfab = x.const_array(mfi);
            auto const& yfab = y.const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                sm += xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

/**
 * \brief Compute dot product of FabArray with itself
 *
 * \param x      FabArray
 * \param xcomp  starting component of x
 * \param ncomp  number of components
 * \param nghost number of ghost cells
 * \param local  If true, MPI communication is skipped.
 */
template <typename FAB, std::enable_if_t<IsBaseFab<FAB>::value,int> FOO = 0>
typename FAB::value_type
Dot (FabArray<FAB> const& x, int xcomp, int ncomp, IntVect const& nghost, bool local = false)
{
    BL_ASSERT(x.nGrowVect().allGE(nghost));

    BL_PROFILE("amrex::Dot()");

    using T = typename FAB::value_type;
    auto sm = T(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& xma = x.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<T>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
        {
            auto t = T(0.0);
            auto const& xfab = xma[box_no];
            for (int n = 0; n < ncomp; ++n) {
                auto v = xfab(i,j,k,xcomp+n);
                t += v*v;
            }
            return t;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& xfab = x.const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                auto v = xfab(i,j,k,xcomp+n);
                sm += v*v;
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

/**
 * \brief Compute dot product of two FabArrays in region that mask is true
 *
 * \param mask   mask
 * \param x      first FabArray
 * \param xcomp  starting component of x
 * \param y      second FabArray
 * \param ycomp  starting component of y
 * \param ncomp  number of components
 * \param nghost number of ghost cells
 * \param local  If true, MPI communication is skipped.
 */
template <typename IFAB, typename FAB,
          std::enable_if_t<IsBaseFab<FAB>::value && IsBaseFab<IFAB>::value,int> FOO = 0>
typename FAB::value_type
Dot (FabArray<IFAB> const& mask, FabArray<FAB> const& x, int xcomp,
     FabArray<FAB> const& y, int ycomp, int ncomp, IntVect const& nghost,
     bool local = false)
{
    BL_ASSERT(x.boxArray() == y.boxArray() && x.boxArray() == mask.boxArray());
    BL_ASSERT(x.DistributionMap() == y.DistributionMap() && x.DistributionMap() == mask.DistributionMap());
    BL_ASSERT(x.nGrowVect().allGE(nghost) && y.nGrowVect().allGE(nghost) &&
              mask.nGrowVect().allGE(nghost));

    BL_PROFILE("amrex::Dot()");

    using T = typename FAB::value_type;
    auto sm = T(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& mma = mask.const_arrays();
        auto const& xma = x.const_arrays();
        auto const& yma = y.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<T>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
        {
            auto t = T(0.0);
            auto m = T(mma[box_no](i,j,k));
            if (m != 0) {
                auto const& xfab = xma[box_no];
                auto const& yfab = yma[box_no];
                for (int n = 0; n < ncomp; ++n) {
                    t += xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
                }
            }
            return t*m;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& mfab = mask.const_array(mfi);
            auto const& xfab = x.const_array(mfi);
            auto const& yfab = y.const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                auto m = T(mfab(i,j,k));
                sm += m * xfab(i,j,k,xcomp+n) * yfab(i,j,k,ycomp+n);
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

/**
 * \brief Compute dot product of FabArray with itself in region that mask is true
 *
 * \param mask   mask
 * \param x      FabArray
 * \param xcomp  starting component of x
 * \param ncomp  number of components
 * \param nghost number of ghost cells
 * \param local  If true, MPI communication is skipped.
 */
template <typename IFAB, typename FAB,
          std::enable_if_t<IsBaseFab<FAB>::value && IsBaseFab<IFAB>::value,int> FOO = 0>
typename FAB::value_type
Dot (FabArray<IFAB> const& mask, FabArray<FAB> const& x, int xcomp, int ncomp,
     IntVect const& nghost, bool local = false)
{
    BL_ASSERT(x.boxArray() == mask.boxArray());
    BL_ASSERT(x.DistributionMap() == mask.DistributionMap());
    BL_ASSERT(x.nGrowVect().allGE(nghost) && mask.nGrowVect().allGE(nghost));

    BL_PROFILE("amrex::Dot()");

    using T = typename FAB::value_type;
    auto sm = T(0.0);
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        auto const& mma = mask.const_arrays();
        auto const& xma = x.const_arrays();
        sm = ParReduce(TypeList<ReduceOpSum>{}, TypeList<T>{}, x, nghost,
        [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept -> GpuTuple<T>
        {
            auto t = T(0.0);
            auto m = T(mma[box_no](i,j,k));
            if (m != 0) {
                auto const& xfab = xma[box_no];
                for (int n = 0; n < ncomp; ++n) {
                    auto v = xfab(i,j,k,xcomp+n);
                    t += v*v;
                }
            }
            return t*m;
        });
    } else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
        for (MFIter mfi(x,true); mfi.isValid(); ++mfi)
        {
            Box const& bx = mfi.growntilebox(nghost);
            auto const& mfab = mask.const_array(mfi);
            auto const& xfab = x.const_array(mfi);
            AMREX_LOOP_4D(bx, ncomp, i, j, k, n,
            {
                auto m = T(mfab(i,j,k));
                auto v = xfab(i,j,k,xcomp+n);
                sm += m*v*v;
            });
        }
    }

    if (!local) {
        ParallelAllReduce::Sum(sm, ParallelContext::CommunicatorSub());
    }

    return sm;
}

//! dst = val
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void setVal (MF& dst, typename MF::value_type val)
{
    dst.setVal(val);
}

//! dst = val in ghost cells.
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void setBndry (MF& dst, typename MF::value_type val, int scomp, int ncomp)
{
    dst.setBndry(val, scomp, ncomp);
}

//! dst *= val
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Scale (MF& dst, typename MF::value_type val, int scomp, int ncomp, int nghost)
{
    dst.mult(val, scomp, ncomp, nghost);
}

//! dst = src
template <class DMF, class SMF,
          std::enable_if_t<IsMultiFabLike_v<DMF> &&
                           IsMultiFabLike_v<SMF>, int> = 0>
void LocalCopy (DMF& dst, SMF const& src, int scomp, int dcomp,
                int ncomp, IntVect const& nghost)
{
    amrex::Copy(dst, src, scomp, dcomp, ncomp, nghost);
}

//! dst += src
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void LocalAdd (MF& dst, MF const& src, int scomp, int dcomp,
                int ncomp, IntVect const& nghost)
{
    amrex::Add(dst, src, scomp, dcomp, ncomp, nghost);
}

//! dst += a * src
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Saxpy (MF& dst, typename MF::value_type a, MF const& src, int scomp, int dcomp,
            int ncomp, IntVect const& nghost)
{
    MF::Saxpy(dst, a, src, scomp, dcomp, ncomp, nghost);
}

//! dst = src + a * dst
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Xpay (MF& dst, typename MF::value_type a, MF const& src, int scomp, int dcomp,
           int ncomp, IntVect const& nghost)
{
    MF::Xpay(dst, a, src, scomp, dcomp, ncomp, nghost);
}

//! dst = a*src_a + b*src_b
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void LinComb (MF& dst,
              typename MF::value_type a, MF const& src_a, int acomp,
              typename MF::value_type b, MF const& src_b, int bcomp,
              int dcomp, int ncomp, IntVect const& nghost)
{
    MF::LinComb(dst, a, src_a, acomp, b, src_b, bcomp, dcomp, ncomp, nghost);
}

//! dst = src w/ MPI communication
template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>, int> = 0>
void ParallelCopy (MF& dst, MF const& src, int scomp, int dcomp, int ncomp,
                   IntVect const& ng_src = IntVect(0),
                   IntVect const& ng_dst = IntVect(0),
                   Periodicity const& period = Periodicity::NonPeriodic())
{
    dst.ParallelCopy(src, scomp, dcomp, ncomp, ng_src, ng_dst, period);
}

template <class MF, std::enable_if_t<IsMultiFabLike_v<MF>, int> = 0>
[[nodiscard]] typename MF::value_type
norminf (MF const& mf, int scomp, int ncomp, IntVect const& nghost,
         bool local = false)
{
    return mf.norminf(scomp, ncomp, nghost, local);
}

//! dst = val
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void setVal (Array<MF,N>& dst, typename MF::value_type val)
{
    for (auto& mf: dst) {
        mf.setVal(val);
    }
}

//! dst = val in ghost cells.
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void setBndry (Array<MF,N>& dst, typename MF::value_type val, int scomp, int ncomp)
{
    for (auto& mf : dst) {
        mf.setBndry(val, scomp, ncomp);
    }
}

//! dst *= val
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Scale (Array<MF,N>& dst, typename MF::value_type val, int scomp, int ncomp,
            int nghost)
{
    for (auto& mf : dst) {
        mf.mult(val, scomp, ncomp, nghost);
    }
}

//! dst = src
template <class DMF, class SMF, std::size_t N,
          std::enable_if_t<IsMultiFabLike_v<DMF> &&
                           IsMultiFabLike_v<SMF>, int> = 0>
void LocalCopy (Array<DMF,N>& dst, Array<SMF,N> const& src, int scomp, int dcomp,
                int ncomp, IntVect const& nghost)
{
    for (std::size_t i = 0; i < N; ++i) {
        amrex::Copy(dst[i], src[i], scomp, dcomp, ncomp, nghost);
    }
}

//! dst += src
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void LocalAdd (Array<MF,N>& dst, Array<MF,N> const& src, int scomp, int dcomp,
               int ncomp, IntVect const& nghost)
{
    for (std::size_t i = 0; i < N; ++i) {
        amrex::Add(dst[i], src[i], scomp, dcomp, ncomp, nghost);
    }
}

//! dst += a * src
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Saxpy (Array<MF,N>& dst, typename MF::value_type a,
            Array<MF,N> const& src, int scomp, int dcomp, int ncomp,
            IntVect const& nghost)
{
    for (std::size_t i = 0; i < N; ++i) {
        MF::Saxpy(dst[i], a, src[i], scomp, dcomp, ncomp, nghost);
    }
}

//! dst = src + a * dst
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void Xpay (Array<MF,N>& dst, typename MF::value_type a,
           Array<MF,N> const& src, int scomp, int dcomp, int ncomp,
           IntVect const& nghost)
{
    for (std::size_t i = 0; i < N; ++i) {
        MF::Xpay(dst[i], a, src[i], scomp, dcomp, ncomp, nghost);
    }
}

//! dst = a*src_a + b*src_b
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>,int> = 0>
void LinComb (Array<MF,N>& dst,
              typename MF::value_type a, Array<MF,N> const& src_a, int acomp,
              typename MF::value_type b, Array<MF,N> const& src_b, int bcomp,
              int dcomp, int ncomp, IntVect const& nghost)
{
    for (std::size_t i = 0; i < N; ++i) {
        MF::LinComb(dst[i], a, src_a[i], acomp, b, src_b[i], bcomp, dcomp, ncomp, nghost);
    }
}

//! dst = src w/ MPI communication
template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>, int> = 0>
void ParallelCopy (Array<MF,N>& dst, Array<MF,N> const& src,
                   int scomp, int dcomp, int ncomp,
                   IntVect const& ng_src = IntVect(0),
                   IntVect const& ng_dst = IntVect(0),
                   Periodicity const& period = Periodicity::NonPeriodic())
{
    for (std::size_t i = 0; i < N; ++i) {
        dst[i].ParallelCopy(src[i], scomp, dcomp, ncomp, ng_src, ng_dst, period);
    }
}

template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF>, int> = 0>
[[nodiscard]] typename MF::value_type
norminf (Array<MF,N> const& mf, int scomp, int ncomp, IntVect const& nghost,
         bool local = false)
{
    auto r = typename MF::value_type(0);
    for (std::size_t i = 0; i < N; ++i) {
        auto tmp = mf[i].norminf(scomp, ncomp, nghost, true);
        r = std::max(r,tmp);
    }
    if (!local) {
        ParallelAllReduce::Max(r, ParallelContext::CommunicatorSub());
    }
    return r;
}

template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF> && (N > 0),
                                                    int> = 0>
[[nodiscard]] int nComp (Array<MF,N> const& mf)
{
    return mf[0].nComp();
}

template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF> && (N > 0),
                                                    int> = 0>
[[nodiscard]] IntVect nGrowVect (Array<MF,N> const& mf)
{
    return mf[0].nGrowVect();
}

template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF> && (N > 0),
                                                    int> = 0>
[[nodiscard]] BoxArray const&
boxArray (Array<MF,N> const& mf)
{
    return mf[0].boxArray();
}

template <class MF, std::size_t N, std::enable_if_t<IsMultiFabLike_v<MF> && (N > 0),
                                                    int> = 0>
[[nodiscard]] DistributionMapping const&
DistributionMap (Array<MF,N> const& mf)
{
    return mf[0].DistributionMap();
}

}

#endif
