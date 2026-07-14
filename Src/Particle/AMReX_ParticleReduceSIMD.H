#ifndef AMREX_PARTICLEREDUCE_SIMD_H_
#define AMREX_PARTICLEREDUCE_SIMD_H_
#include <AMReX_Config.H>

#include <AMReX_ParticleReduce.H>  // provides amrex::SIMDindex via AMReX_Gpu.H
#include <AMReX_ReduceSIMD.H>
#include <AMReX_SIMD.H>
#include <AMReX_Tuple.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Vector.H>

#include <type_traits>
#include <utility>

namespace amrex {

/// \cond DOXYGEN_IGNORE
namespace particle_detail {

    //! Wraps a callable that returns a single value into one that returns a 1-element tuple
    template <typename F>
    struct SIMDTupleWrapper
    {
        F m_f;

        template <typename PTD, typename SI>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        auto operator() (PTD const& ptd, SI const si) const noexcept
        {
            return amrex::makeTuple(m_f(ptd, si));
        }
    };
}
/// \endcond

/**
 * \brief A general SIMD-vectorized reduction method for the particles in a
 * ParticleContainer. This version operates from the specified lev_min to lev_max.
 *
 * This is the SIMD counterpart of amrex::ParticleReduce. When AMReX is compiled
 * with SIMD support (\c AMReX_SIMD=ON), the reduction runs vectorized over
 * \c WIDTH particles at a time (with a scalar remainder loop). The same user
 * code compiles and runs unchanged in all other configurations:
 *   - \c AMReX_SIMD=OFF (CPU): a scalar loop, identical in semantics to
 *     amrex::ParticleReduce.
 *   - GPU builds: the regular device reduction, like amrex::ParticleReduce.
 *
 * The callable \c f is invoked as \c f(ptd,si) with a \c ConstParticleTileData
 * and an amrex::SIMDindex. Write it generically over the index and load
 * components with amrex::simd::load_1d so that a single source supports the
 * vectorized main loop, the scalar remainder and GPUs:
 *
 *    \code{.cpp}
 *       using ConstPTDType = typename PC::ConstPTDType;
 *       amrex::ReduceOps<amrex::ReduceOpSum, amrex::ReduceOpMin> reduce_ops;
 *       auto r = amrex::ParticleReduceSIMD<amrex::ReduceData<amrex::ParticleReal, amrex::ParticleReal>>(
 *                    pc, lev_min, lev_max,
 *                    [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, auto si) noexcept
 *                {
 *                    auto const a = amrex::simd::load_1d(ptd.rdata(1), si);
 *                    auto const w = amrex::simd::load_1d(ptd.rdata(2), si);
 *                    return amrex::makeTuple(a * w, a);
 *                }, reduce_ops);
 *    \endcode
 *
 * For CUDA builds, note that extended device lambdas have compiler-version
 * dependent restrictions on generic (auto) parameters; a small functor with a
 * templated \c operator() is the most portable form of the kernel.
 *
 * Only \c ReduceOpSum, \c ReduceOpMin and \c ReduceOpMax are supported
 * (compile-time enforced). Note that there is no MPI reduction performed at
 * the end of this operation, just like for amrex::ParticleReduce.
 *
 * \note With SIMD enabled, sum reductions reassociate floating-point additions
 *       (per-lane partial sums plus a final horizontal fold), so results can
 *       differ from the scalar evaluation order by rounding. In particular,
 *       \c amrex::system::regtest_reduction is not honored by the SIMD path.
 *
 * \tparam RD the ReduceData type
 * \tparam WIDTH SIMD width in elements; defaults to the native SIMD register
 *         width for amrex::ParticleReal (1 without SIMD support)
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 * \tparam ReduceOps_t the ReduceOps type
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a callable that operates on WIDTH particles at a time, see above
 * \param reduce_ops the reduction operations, one per component of RD
 */
template <class RD, int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F, class ReduceOps_t>
requires (IsParticleContainer<PC>::value)
typename RD::Type
ParticleReduceSIMD (PC const& pc, int lev_min, int lev_max, F const& f, ReduceOps_t& reduce_ops)
{
    RD reduce_data(reduce_ops);
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        const auto& plev = pc.GetParticles(lev);
        Vector<std::pair<int, int> > grid_tile_ids;
        Vector<const typename PC::ParticleTileType*> ptile_ptrs;
        for (auto& kv : plev)
        {
            grid_tile_ids.push_back(kv.first);
            ptile_ptrs.push_back(&(kv.second));
        }

        // We avoid ParIter here to allow reductions after grids change but prior to calling Redistribute()
#if !defined(AMREX_USE_GPU) && defined(AMREX_USE_OMP)
#pragma omp parallel for
#endif
        for (int pmap_it = 0; pmap_it < std::ssize(ptile_ptrs); ++pmap_it)
        {
            const auto& tile = plev.at(grid_tile_ids[pmap_it]);
            const auto np = tile.numParticles();
            const auto& ptd = tile.getConstParticleTileData();
#ifdef AMREX_USE_GPU
            reduce_ops.eval(np, reduce_data,
                            [=] AMREX_GPU_DEVICE (const int i) noexcept
                            {
                                return f(ptd, SIMDindex<1, int>{i});
                            });
#else
#ifdef AMREX_USE_SIMD
            constexpr int EffWidth = WIDTH;
#else
            constexpr int EffWidth = 1;  // explicit WIDTH > 1 still compiles, scalar
#endif
            evalReduceSIMD<EffWidth>(np, reduce_data,
                            [=] (auto si) noexcept
                            {
                                return f(ptd, si);
                            }, reduce_ops);
#endif
        }
    }
    return reduce_data.value(reduce_ops);
}

/**
 * \brief A general SIMD-vectorized reduction method for the particles in a
 * ParticleContainer. This version operates only on the specified level.
 *
 * See the lev_min/lev_max overload of amrex::ParticleReduceSIMD for the
 * kernel contract, an example, and the supported reduce operations.
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a callable that operates on WIDTH particles at a time
 * \param reduce_ops the reduction operations, one per component of RD
 */
template <class RD, int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F, class ReduceOps_t>
requires (IsParticleContainer<PC>::value)
typename RD::Type
ParticleReduceSIMD (PC const& pc, int lev, F const& f, ReduceOps_t& reduce_ops)
{
    return ParticleReduceSIMD<RD, WIDTH>(pc, lev, lev, f, reduce_ops);
}

/**
 * \brief A general SIMD-vectorized reduction method for the particles in a
 * ParticleContainer. This version operates over all particles on all levels.
 *
 * See the lev_min/lev_max overload of amrex::ParticleReduceSIMD for the
 * kernel contract, an example, and the supported reduce operations.
 *
 * \param pc the ParticleContainer to operate on
 * \param f a callable that operates on WIDTH particles at a time
 * \param reduce_ops the reduction operations, one per component of RD
 */
template <class RD, int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F, class ReduceOps_t>
requires (IsParticleContainer<PC>::value)
typename RD::Type
ParticleReduceSIMD (PC const& pc, F const& f, ReduceOps_t& reduce_ops)
{
    return ParticleReduceSIMD<RD, WIDTH>(pc, 0, pc.finestLevel(), f, reduce_ops);
}

/**
 * \brief A SIMD-vectorized sum reduction over the particles in a ParticleContainer.
 * This version operates from the specified lev_min to lev_max.
 *
 * This is the SIMD counterpart of the \c (ptd,i) form of amrex::ReduceSum. The
 * callable is invoked as \c f(ptd,si) with an amrex::SIMDindex and must return
 * one value (a SIMD register of \c WIDTH values in the vectorized main loop, a
 * scalar in the remainder loop and on GPUs). Load particle components with
 * amrex::simd::load_1d so a single source supports all configurations:
 *
 *    \code{.cpp}
 *       using ConstPTDType = typename PC::ConstPTDType;
 *       auto sm = amrex::ReduceSumSIMD(pc, lev_min, lev_max,
 *                 [=] AMREX_GPU_DEVICE (const ConstPTDType& ptd, auto si) noexcept
 *                 {
 *                     return amrex::simd::load_1d(ptd.rdata(0), si);
 *                 });
 *    \endcode
 *
 * Note that there is no MPI reduction performed at the end of this operation.
 * On GPU builds this always runs on the device (like amrex::ParticleReduce and
 * unlike amrex::ReduceSum, which checks the launch region).
 *
 * \note With SIMD enabled, sum reductions reassociate floating-point additions,
 *       so results can differ from the scalar evaluation order by rounding.
 *
 * \tparam WIDTH SIMD width in elements; defaults to the native SIMD register
 *         width for amrex::ParticleReal (1 without SIMD support)
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a callable that operates on WIDTH particles at a time, see above
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceSumSIMD (PC const& pc, int lev_min, int lev_max, F const& f)
{
    using ConstPTDType = typename PC::ConstPTDType;
    using T = std::decay_t<std::invoke_result_t<F, ConstPTDType const&, SIMDindex<1, int>>>;
    ReduceOps<ReduceOpSum> reduce_ops;
    auto r = ParticleReduceSIMD<ReduceData<T>, WIDTH>(
        pc, lev_min, lev_max, particle_detail::SIMDTupleWrapper<F>{f}, reduce_ops);
    return amrex::get<0>(r);
}

/**
 * \brief A SIMD-vectorized sum reduction over the particles in a ParticleContainer.
 * This version operates only on the specified level.
 *
 * See the lev_min/lev_max overload of amrex::ReduceSumSIMD for the kernel
 * contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceSumSIMD (PC const& pc, int lev, F const& f)
{
    return ReduceSumSIMD<WIDTH>(pc, lev, lev, f);
}

/**
 * \brief A SIMD-vectorized sum reduction over the particles in a ParticleContainer.
 * This version operates over all particles on all levels.
 *
 * See the lev_min/lev_max overload of amrex::ReduceSumSIMD for the kernel
 * contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceSumSIMD (PC const& pc, F const& f)
{
    return ReduceSumSIMD<WIDTH>(pc, 0, pc.finestLevel(), f);
}

/**
 * \brief A SIMD-vectorized min reduction over the particles in a ParticleContainer.
 * This version operates from the specified lev_min to lev_max.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example; this
 * version uses "min" as the reduction operation (exact, no floating-point
 * reassociation concerns).
 *
 * \tparam WIDTH SIMD width in elements; defaults to the native SIMD register
 *         width for amrex::ParticleReal (1 without SIMD support)
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMinSIMD (PC const& pc, int lev_min, int lev_max, F const& f)
{
    using ConstPTDType = typename PC::ConstPTDType;
    using T = std::decay_t<std::invoke_result_t<F, ConstPTDType const&, SIMDindex<1, int>>>;
    ReduceOps<ReduceOpMin> reduce_ops;
    auto r = ParticleReduceSIMD<ReduceData<T>, WIDTH>(
        pc, lev_min, lev_max, particle_detail::SIMDTupleWrapper<F>{f}, reduce_ops);
    return amrex::get<0>(r);
}

/**
 * \brief A SIMD-vectorized min reduction over the particles in a ParticleContainer.
 * This version operates only on the specified level.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMinSIMD (PC const& pc, int lev, F const& f)
{
    return ReduceMinSIMD<WIDTH>(pc, lev, lev, f);
}

/**
 * \brief A SIMD-vectorized min reduction over the particles in a ParticleContainer.
 * This version operates over all particles on all levels.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMinSIMD (PC const& pc, F const& f)
{
    return ReduceMinSIMD<WIDTH>(pc, 0, pc.finestLevel(), f);
}

/**
 * \brief A SIMD-vectorized max reduction over the particles in a ParticleContainer.
 * This version operates from the specified lev_min to lev_max.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example; this
 * version uses "max" as the reduction operation (exact, no floating-point
 * reassociation concerns).
 *
 * \tparam WIDTH SIMD width in elements; defaults to the native SIMD register
 *         width for amrex::ParticleReal (1 without SIMD support)
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMaxSIMD (PC const& pc, int lev_min, int lev_max, F const& f)
{
    using ConstPTDType = typename PC::ConstPTDType;
    using T = std::decay_t<std::invoke_result_t<F, ConstPTDType const&, SIMDindex<1, int>>>;
    ReduceOps<ReduceOpMax> reduce_ops;
    auto r = ParticleReduceSIMD<ReduceData<T>, WIDTH>(
        pc, lev_min, lev_max, particle_detail::SIMDTupleWrapper<F>{f}, reduce_ops);
    return amrex::get<0>(r);
}

/**
 * \brief A SIMD-vectorized max reduction over the particles in a ParticleContainer.
 * This version operates only on the specified level.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMaxSIMD (PC const& pc, int lev, F const& f)
{
    return ReduceMaxSIMD<WIDTH>(pc, lev, lev, f);
}

/**
 * \brief A SIMD-vectorized max reduction over the particles in a ParticleContainer.
 * This version operates over all particles on all levels.
 *
 * See amrex::ReduceSumSIMD for the kernel contract and an example.
 *
 * \param pc the ParticleContainer to operate on
 * \param f a callable that operates on WIDTH particles at a time
 */
template <int WIDTH = static_cast<int>(simd::native_simd_size_particlereal),
          class PC, class F>
requires (IsParticleContainer<PC>::value)
auto
ReduceMaxSIMD (PC const& pc, F const& f)
{
    return ReduceMaxSIMD<WIDTH>(pc, 0, pc.finestLevel(), f);
}

}

#endif
