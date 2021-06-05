#ifndef AMREX_PARTICLEREDUCE_H_
#define AMREX_PARTICLEREDUCE_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_GpuUtility.H>
#include <AMReX_TypeTraits.H>

#include <limits>

namespace amrex
{

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version uses "Sum" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceSum (PC const& pc, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceSum(pc, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version uses "Sum" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceSum (PC const& pc, int lev, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceSum(pc, lev, lev, std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version uses "Sum" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceSum (PC const& pc, int lev_min, int lev_max, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    using value_type = decltype(f(typename PC::SuperParticleType()));
    using ParIter = typename PC::ParConstIterType;
    value_type sm = 0;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<value_type> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto ptd = tile.getConstParticleTileData();
                reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple {return {f(ptd.getSuperParticle(i))};});
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        sm = amrex::get<0>(hv);
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(+:sm)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                for (int i = 0; i < np; ++i)
                    sm += f(ptd.getSuperParticle(i));
            }
        }
    }

    return sm;
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version uses "Max" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMax (PC const& pc, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceMax(pc, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version uses "Mas" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMax (PC const& pc, int lev, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceMax(pc, lev, lev, std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version uses "Max" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMax (PC const& pc, int lev_min, int lev_max, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    using value_type = decltype(f(typename PC::SuperParticleType()));
    using ParIter = typename PC::ParConstIterType;
    constexpr value_type value_lowest = std::numeric_limits<value_type>::lowest();
    value_type r = value_lowest;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpMax> reduce_op;
        ReduceData<value_type> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple {return {f(ptd.getSuperParticle(i))};});
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        r = amrex::get<0>(hv);
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(max:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                for (int i = 0; i < np; ++i)
                    r = std::max(r, f(ptd.getSuperParticle(i)));
            }
        }
    }

    return r;
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version uses "Min" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMin (PC const& pc, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceMin(pc, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version uses "Min" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMin (PC const& pc, int lev, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    return ReduceMin(pc, lev, lev, std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version uses "Min" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
auto
ReduceMin (PC const& pc, int lev_min, int lev_max, F&& f) -> decltype(f(typename PC::SuperParticleType()))
{
    using value_type = decltype(f(typename PC::SuperParticleType()));
    using ParIter = typename PC::ParConstIterType;
    constexpr value_type value_max = std::numeric_limits<value_type>::max();
    value_type r = value_max;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpMin> reduce_op;
        ReduceData<value_type> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple {return {f(ptd.getSuperParticle(i))};});
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        r = amrex::get<0>(hv);
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(min:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                for (int i = 0; i < np; ++i)
                    r = std::min(r, f(ptd.getSuperParticle(i)));
            }
        }
    }

    return r;
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version uses "LogicalAnd" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalAnd (PC const& pc, F&& f)
{
    return ReduceLogicalAnd(pc, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version uses "LogicalAnd" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalAnd (PC const& pc, int lev, F&& f)
{
    return ReduceLogicalAnd(pc, lev, lev, std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version uses "LogicalAnd" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalAnd (PC const& pc, int lev_min, int lev_max, F&& f)
{
    using ParIter = typename PC::ParConstIterType;
    int r = true;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpLogicalAnd> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                reduce_op.eval(np, reduce_data,
                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple {return {f(ptd.getSuperParticle(i))};});
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        r = amrex::get<0>(hv);
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(&&:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                for (int i = 0; i < np; ++i)
                    r = r && f(ptd.getSuperParticle(i));
            }
        }
    }

    return r;
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version uses "LogicalOr" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalOr (PC const& pc, F&& f)
{
    return ReduceLogicalOr(pc, 0, pc.finestLevel(), std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version uses "LogicalOr" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalOr (PC const& pc, int lev, F&& f)
{
    return ReduceLogicalOr(pc, lev, lev, std::forward<F>(f));
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version uses "LogicalOr" as the reduction operation. The quantity reduced over is an arbitrary function
 * of a "superparticle", which contains all the data in the particle type, whether it is stored in AoS or
 * SoA form.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns the value to be reduced over all particles.
 *
 */
template <class PC, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
bool
ReduceLogicalOr (PC const& pc, int lev_min, int lev_max, F&& f)
{
    using ParIter = typename PC::ParConstIterType;
    int r = false;

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        ReduceOps<ReduceOpLogicalOr> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                reduce_op.eval(np, reduce_data,
                                [=] AMREX_GPU_DEVICE (const int i) -> ReduceTuple {return {f(ptd.getSuperParticle(i))};});
            }
        }

        ReduceTuple hv = reduce_data.value(reduce_op);
        r = amrex::get<0>(hv);
    }
    else
#endif
    {
        for (int lev = lev_min; lev <= lev_max; ++lev)
        {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (!system::regtest_reduction) reduction(||:r)
#endif
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                                const auto ptd = tile.getConstParticleTileData();
                for (int i = 0; i < np; ++i)
                    r = r || f(ptd.getSuperParticle(i));

            }
        }
    }

    return r;
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates over all particles on all levels.
 *
 * This version can operate on a GpuTuple worth of data at once.
 * It also takes an arbitrary tuple of reduction operators.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * Unlike the other reduction functions in this file, this version does not respect the
 * Gpu::launchRegion flag. If AMReX is built with GPU support, this reduction will always be
 * done on the device.
 *
 * \tparam RD an amrex::ReduceData type
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 * \tparam ReduceOps a ReduceOps type
 *
 * \param pc the ParticleContainer to operate on
 * \param f a function that takes a "superparticle" and returns a tuple of values to be reduced over all particles.
 * \param reduce_ops specifies the reduction operations for each tuple element
 *
 * Example usage:
 *    amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
 *    auto r = amrex::ParticleReduce<ReduceData<amrex::Real, amrex::Real,int>> (
 *                 pc, [=] AMREX_GPU_DEVICE (const PType& p) noexcept
 *                               -> amrex::GpuTuple<amrex::Real,amrex::Real,int>
 *             {
 *                 const amrex::Real a = p.rdata(1);
 *                 const amrex::Real b = p.rdata(2);
 *                 const int c = p.idata(1);
 *                 return {a, b, c};
 *             }, reduce_ops);
 *
 */
template <class RD, class PC, class F, class ReduceOps,
          std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
typename RD::Type
ParticleReduce (PC const& pc, F&& f, ReduceOps& reduce_ops)
{
    return ParticleReduce<RD>(pc, 0, pc.finestLevel(), std::forward<F>(f), reduce_ops);
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates only on the specified level.
 *
 * This version can operate on a GpuTuple worth of data at once.
 * It also takes an arbitrary tuple of reduction operators.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * Unlike the other reduction functions in this file, this version does not respect the
 * Gpu::launchRegion flag. If AMReX is built with GPU support, this reduction will always be
 * done on the device.
 *
 * \tparam RD an amrex::ReduceData type
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 * \tparam ReduceOps a ReduceOps type
 *
 * \param pc the ParticleContainer to operate on
 * \param lev the level to operate on
 * \param f a function that takes a "superparticle" and returns a tuple of values to be reduced over all particles.
 * \param reduce_ops specifies the reduction operations for each tuple element
 *
 * Example usage:
 *    amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
 *    auto r = amrex::ParticleReduce<ReduceData<amrex::Real, amrex::Real,int>> (
 *                 pc, 0, [=] AMREX_GPU_DEVICE (const PType& p) noexcept
 *                               -> amrex::GpuTuple<amrex::Real,amrex::Real,int>
 *             {
 *                 const amrex::Real a = p.rdata(1);
 *                 const amrex::Real b = p.rdata(2);
 *                 const int c = p.idata(1);
 *                 return {a, b, c};
 *             }, reduce_ops);
 *
 */
template <class RD, class PC, class F, class ReduceOps,
          std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
typename RD::Type
ParticleReduce (PC const& pc, int lev, F&& f, ReduceOps& reduce_ops)
{
    return ParticleReduce<RD>(pc, lev, lev, std::forward<F>(f), reduce_ops);
}

/**
 * \brief A general reduction method for the particles in a ParticleContainer that can run on either CPUs or GPUs.
 * This version operates from the specified lev_min to lev_max.
 *
 * This version can operate on a GpuTuple worth of data at once.
 * It also takes an arbitrary tuple of reduction operators.
 *
 * Note that there is no MPI reduction performed at the end of this operation. Users should manually
 * call the MPI reduction operations described in ParallelDescriptor if they want that behavior.
 *
 * Unlike the other reduction functions in this file, this version does not respect the
 * Gpu::launchRegion flag. If AMReX is built with GPU support, this reduction will always be
 * done on the device.
 *
 * \tparam RD an amrex::ReduceData type
 * \tparam PC the ParticleContainer type
 * \tparam F a function object
 * \tparam ReduceOps a ReduceOps type
 *
 * \param pc the ParticleContainer to operate on
 * \param lev_min the minimum level to include
 * \param lev_max the maximum level to include
 * \param f a function that takes a "superparticle" and returns a tuple of values to be reduced over all particles.
 * \param reduce_ops specifies the reduction operations for each tuple element
 *
 * Example usage:
 *    amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
 *    auto r = amrex::ParticleReduce<ReduceData<amrex::Real, amrex::Real,int>> (
 *                 pc, 0, 1, [=] AMREX_GPU_DEVICE (const PType& p) noexcept
 *                               -> amrex::GpuTuple<amrex::Real,amrex::Real,int>
 *             {
 *                 const amrex::Real a = p.rdata(1);
 *                 const amrex::Real b = p.rdata(2);
 *                 const int c = p.idata(1);
 *                 return {a, b, c};
 *             }, reduce_ops);
 *
 */
template <class RD, class PC, class F, class ReduceOps,
          std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
typename RD::Type
ParticleReduce (PC const& pc, int lev_min, int lev_max, F&& f, ReduceOps& reduce_ops)
{
    RD reduce_data(reduce_ops);
    using ParIter = typename PC::ParConstIterType;
#if !defined(AMREX_USE_GPU) && defined(AMREX_USE_OMP)
#pragma omp parallel
#endif
    for (int lev = lev_min; lev <= lev_max; ++lev) {
        for (ParIter pti(pc, lev); pti.isValid(); ++pti) {
            const auto& tile = pti.GetParticleTile();
            const auto np = tile.numParticles();
            const auto ptd = tile.getConstParticleTileData();
            reduce_ops.eval(np, reduce_data,
                            [=] AMREX_GPU_DEVICE (const int i) noexcept
                            {
                                return f(ptd.getSuperParticle(i));
                            });
        }
    }
    return reduce_data.value(reduce_ops);
}
}
#endif
