#ifndef AMREX_PARTICLEMESH_H_
#define AMREX_PARTICLEMESH_H_
#include <AMReX_Config.H>

#include <AMReX_TypeTraits.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParticleUtil.H>
#include <type_traits>

namespace amrex {

namespace particle_detail {

template <typename F, typename T, typename T_ParticleType, template<class, int, int> class PTDType, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const PTDType<T_ParticleType, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
{
    using PTDTypeT = std::remove_const_t<std::remove_reference_t<decltype(p)>>;
    if constexpr ( ! T_ParticleType::is_soa_particle &&
                   IsCallable<F, typename PTDTypeT::ParticleRefType, decltype(fabarr), decltype(plo), decltype(dxi)>::value) {
        return f(p.m_aos[i], fabarr, plo, dxi);
    } else if constexpr ( ! T_ParticleType::is_soa_particle &&
                          IsCallable<F, typename PTDTypeT::ParticleRefType, decltype(fabarr)>::value) {
        return f(p.m_aos[i], fabarr);
    } else if constexpr (IsCallable<F, decltype(p.getSuperParticle(i)), decltype(fabarr), decltype(plo), decltype(dxi)>::value) {
        return f(p.getSuperParticle(i), fabarr, plo, dxi);
    } else if constexpr (IsCallable<F, decltype(p.getSuperParticle(i)), decltype(fabarr)>::value) {
        return f(p.getSuperParticle(i), fabarr);
    } else if constexpr (IsCallable<F, decltype(p), int, decltype(fabarr), decltype(plo), decltype(dxi)>::value) {
        return f(p, i, fabarr, plo, dxi);
    } else {
        return f(p, i, fabarr);
    }
}
}

template <class PC, class MF, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void
ParticleToMesh (PC const& pc, MF& mf, int lev, F const& f, bool zero_out_input=true)
{
    BL_PROFILE("amrex::ParticleToMesh");

    if (zero_out_input) { mf.setVal(0.0); }

    MF* mf_pointer;

    if (pc.OnSameGrids(lev, mf) && zero_out_input)
    {
        mf_pointer = &mf;
    } else {
        mf_pointer = new MF(pc.ParticleBoxArray(lev),
                            pc.ParticleDistributionMap(lev),
                            mf.nComp(), mf.nGrowVect());
        mf_pointer->setVal(0.0);
    }

    const auto plo = pc.Geom(lev).ProbLoArray();
    const auto dxi = pc.Geom(lev).InvCellSizeArray();

    using ParIter = typename PC::ParConstIterType;
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion())
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto& tile = pti.GetParticleTile();
            const auto np = tile.numParticles();
            const auto& ptd = tile.getConstParticleTileData();

            auto& fab = (*mf_pointer)[pti];
            auto fabarr = fab.array();

            AMREX_FOR_1D( np, i,
            {
                particle_detail::call_f(f, ptd, i, fabarr, plo, dxi);
            });
        }
    }
    else
#endif
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
        {
            typename MF::FABType::value_type local_fab;
            for(ParIter pti(pc, lev); pti.isValid(); ++pti)
            {
                const auto& tile = pti.GetParticleTile();
                const auto np = tile.numParticles();
                const auto& ptd = tile.getConstParticleTileData();

                auto& fab = (*mf_pointer)[pti];

                Box tile_box = pti.tilebox();
                tile_box.grow(mf_pointer->nGrowVect());
                local_fab.resize(tile_box,mf_pointer->nComp());
                local_fab.template setVal<RunOn::Host>(0.0);
                auto fabarr = local_fab.array();

                AMREX_FOR_1D( np, i,
                {
                    particle_detail::call_f(f, ptd, i, fabarr, plo, dxi);
                });

                fab.template atomicAdd<RunOn::Host>(local_fab, tile_box, tile_box,
                                                    0, 0, mf_pointer->nComp());
            }
        }
    }

    if (mf_pointer != &mf)
    {
        mf.ParallelAdd(*mf_pointer, 0, 0, mf_pointer->nComp(),
                       mf_pointer->nGrowVect(), IntVect(0), pc.Geom(lev).periodicity());
        delete mf_pointer;
    } else {
        mf_pointer->SumBoundary(pc.Geom(lev).periodicity());
    }
}

template <class PC, class MF, class F, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
void
MeshToParticle (PC& pc, MF const& mf, int lev, F const& f)
{
    BL_PROFILE("amrex::MeshToParticle");

    MF* mf_pointer = pc.OnSameGrids(lev, mf) ?
        const_cast<MF*>(&mf) : new MF(pc.ParticleBoxArray(lev),
                                      pc.ParticleDistributionMap(lev),
                                      mf.nComp(), mf.nGrowVect());

    if (mf_pointer != &mf) {mf_pointer->ParallelCopy(mf,0,0,mf.nComp(),mf.nGrowVect(),mf.nGrowVect()); }

    const auto plo = pc.Geom(lev).ProbLoArray();
    const auto dxi = pc.Geom(lev).InvCellSizeArray();

    using ParIter = typename PC::ParIterType;
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
    for(ParIter pti(pc, lev); pti.isValid(); ++pti)
    {
        auto& tile = pti.GetParticleTile();
        const auto np = tile.numParticles();
        const auto& ptd = tile.getParticleTileData();

        const auto& fab = (*mf_pointer)[pti];
        auto fabarr = fab.array();

        AMREX_FOR_1D( np, i,
        {
            particle_detail::call_f(f, ptd, i, fabarr, plo, dxi);
        });
    }

    if (mf_pointer != &mf) { delete mf_pointer; }
}

}
#endif
