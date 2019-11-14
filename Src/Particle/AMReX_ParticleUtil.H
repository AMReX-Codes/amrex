#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParticleTile.H>

#include <limits>

namespace amrex
{

/**
 * /brief A general single particle copying routine that can run on the GPU. 
 *
 * \tparam NSR number of extra reals in the particle struct
 * \tparam NSI number of extra ints in the particle struct
 * \tparam NAR number of reals in the struct-of-arrays
 * \tparam NAI number of ints in the struct-of-arrays
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_i the index in the source to read from
 * \param dst_i the index in the destination to write to
 *
 */
template <int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_INLINE
void copyParticle(const      ParticleTileData<NSR, NSI, NAR, NAI>& dst, 
                  const ConstParticleTileData<NSR, NSI, NAR, NAI>& src, 
                  int src_i, int dst_i) noexcept
{
    AMREX_ASSERT(dst.m_num_runtime_real == src.m_num_runtime_real);
    AMREX_ASSERT(dst.m_num_runtime_int  == src.m_num_runtime_int );

    dst.m_aos[dst_i] = src.m_aos[src_i];
    for (int j = 0; j < NAR; ++j)
        dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
    for (int j = 0; j < dst.m_num_runtime_real; ++j)
        dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
    for (int j = 0; j < NAI; ++j)
        dst.m_idata[j][dst_i] = src.m_idata[j][src_i];
    for (int j = 0; j < dst.m_num_runtime_int; ++j)
        dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
}

/**
 * /brief Apply the function f to all the particles in src, writing the
 * result to dst. This version does all the particles in src.
 * 
 * \tparam PTile the particle tile type
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */    
template <typename PTile, typename F>
void transformParticles(PTile& dst, const PTile& src, F f) noexcept
{
    auto np = src.numParticles();
    transformParticles(dst, src, 0, 0, np, f);
}

/**
 * /brief Apply the function f to particles in src, writing the
 * result to dst. This version applies the function to n particles
 * starting at index src_start, writing the result starting at src_start.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam N the size type, e.g. long
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param f the function that will be applied to each particle
 *
 */
template <typename PTile, typename IndexType, typename N, typename F>
void transformParticles(PTile& dst, const PTile& src, IndexType src_start, IndexType dst_start, N n, F f) noexcept
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_FOR_1D( n, i, { f(dst_data, src_data, src_start+i, dst_start+i); });
}

/**
 * /brief Conditionally copy particles from src to dst based on the value of mask.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 *
 */
template <typename PTile, typename IndexType>
void filterParticles(PTile& dst, const PTile& src, const IndexType* mask) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets);
    auto p_offsets = offsets.dataPtr();
    
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_FOR_1D( np, i,
    {
        if (mask[i]) copyParticle(dst_data, src_data, i, p_offsets[i]);
    });
}

/**
 * /brief Conditionally copy particles from src to dst based on a predicate.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam Pred a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 *
 */    
template <typename PTile, typename IndexType, typename Pred>
void filterParticles(PTile& dst, const PTile& src, Pred p) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();
    
    AMREX_FOR_1D(np, i, {p_mask[i] = p(src_data, i); });

    filterPartcles(dst, src, mask.dataPtr());
}

/**
 * /brief Conditionally copy particles from src to dst based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename PTile, typename IndexType, typename F>
void filterAndTransformParticles(PTile& dst, const PTile& src, const IndexType* mask, F f) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets);
    auto p_offsets = offsets.dataPtr();
    
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_FOR_1D( np, i,
    {
        if (mask[i]) f(dst_data, src_data, i, p_offsets[i]);
    });
}

/**
 * /brief Conditionally copy particles from src to dst based on a predicate.
 * A transformation will also be applied to the particles on copy.
 * 
 * \tparam PTile the particle tile type
 * \tparam IndexType the index type, e.g. unsigned int
 * \tparam Pred a function object
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */    
template <typename PTile, typename IndexType, typename Pred, typename F>
void filterAndTransformParticles(PTile& dst, const PTile& src, Pred p, F f) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<IndexType> mask(np);

    auto p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();
    
    AMREX_FOR_1D(np, i, {p_mask[i] = p(src_data, i); });

    filterAndTransformPartcles(dst, src, mask.dataPtr(), f);
}
    
/**
 * /brief Gather particles copies particles into contiguous order from an 
 * arbitrary order. Specifically, the particle at the index inds[i] in src
 * will be copied to the index i in dst.
 * 
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. long
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename IndexType>
void gatherParticles (PTile& dst, const PTile& src, N np, const IndexType* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();
          
          AMREX_FOR_1D( np, i, { copyParticle(dst_data, src_data, inds[i], i); });
}

/**
 * /brief Scatter particles copies particles from contiguous order into an
 * arbitrary order. Specifically, the particle at the index i in src
 * will be copied to the index inds[i] in dst.
 * 
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. long
 * \tparam IndexType the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename IndexType>
void scatterParticles (PTile& dst, const PTile& src, N np, const IndexType* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();
          
          AMREX_FOR_1D( np, i, { copyParticle(dst_data, src_data, i, inds[i]); });
}
    
AMREX_GPU_HOST_DEVICE
int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
                  const IntVect& a_tile_size, Box& tbx);

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getParticleCell (P const& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    IntVect iv = IntVect(
        AMREX_D_DECL(floor((p.pos(0)-plo[0])*dxi[0]),
                     floor((p.pos(1)-plo[1])*dxi[1]),
                     floor((p.pos(2)-plo[2])*dxi[2])));
    iv += domain.smallEnd();
    return iv;
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getParticleGrid (P const& p, amrex::Array4<int> const& mask,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                     const Box& domain) noexcept
{
    if (p.id() < 0) return -1;    
    IntVect iv = getParticleCell(p, plo, dxi, domain);
    return mask(iv);
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void enforcePeriodic (P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& phi,
                      amrex::GpuArray<int,AMREX_SPACEDIM> const& is_per) noexcept
{
    // In rare cases, degenerate cases can be found, filter out with a tolerance
#if !defined(__CUDACC__) || (__CUDACC_VER_MAJOR__ != 9) || (__CUDACC_VER_MINOR__ != 2)
    static constexpr Real eps = std::numeric_limits<Real>::epsilon();
#else
    static constexpr Real eps = DBL_EPSILON;
#endif

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (not is_per[idim]) continue;
        if (p.pos(idim) >= phi[idim]) { 
            while (p.pos(idim) >= phi[idim]) { 
                p.pos(idim) -= (phi[idim] - plo[idim]);
            }
            if (p.pos(idim) < plo[idim]) p.pos(idim) = plo[idim]; // clamp to avoid precision issues;
        } 
        else if (p.pos(idim) < plo[idim]) {
            while (p.pos(idim) < plo[idim]) {
                p.pos(idim) += (phi[idim] - plo[idim]); 
            }
            if (p.pos(idim) == phi[idim]) p.pos(idim) = plo[idim]; // clamp to avoid precision issues;
            if (p.pos(idim) > phi[idim]) p.pos(idim) = phi[idim]-eps; // clamp to avoid precision issues;
        }
    }
}

template <class PC>
bool
numParticlesOutOfRange (PC const& pc, int nGrow)
{
    return numParticlesOutOfRange(pc, 0, pc.finestLevel(), nGrow);
}

template <class PC>
bool
numParticlesOutOfRange (PC const& pc, int lev_min, int lev_max, int nGrow)
{
    BL_PROFILE("numParticlesOutOfRange()");

    using ParIter = typename PC::ParConstIterType;
    using ParticleType = typename PC::ParticleType;
    int num_wrong = 0;
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
        const auto plo = pc.Geom(lev).ProbLoArray();
        const auto dxi = pc.Geom(lev).InvCellSizeArray();
        const Box domain = pc.Geom(lev).Domain();
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion() and !system::regtest_reduction) reduction(+:num_wrong)
#endif
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            const auto& tile = pti.GetParticleTile();
            const auto np = tile.numParticles();
            const auto& aos = tile.GetArrayOfStructs();
            const auto pstruct = aos().dataPtr();
            
            Box box = pti.tilebox();
            box.grow(nGrow);

            ReduceOps<ReduceOpSum> reduce_op;
            ReduceData<int> reduce_data(reduce_op);
            using ReduceTuple = typename decltype(reduce_data)::Type;

            reduce_op.eval(np, reduce_data,
            [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
            {
                const ParticleType& p = pstruct[i];
                if ((p.id() < 0)) return false;
                IntVect iv = IntVect(
                        AMREX_D_DECL(floor((p.pos(0)-plo[0])*dxi[0]),
                                     floor((p.pos(1)-plo[1])*dxi[1]),
                                     floor((p.pos(2)-plo[2])*dxi[2])));
                iv += domain.smallEnd();
                return !box.contains(iv);
            });
            int hv = amrex::get<0>(reduce_data.value());
            num_wrong += hv;
        }
    }
    ParallelDescriptor::ReduceIntSum(num_wrong);

    return num_wrong;
}

IntVect computeRefFac (const ParGDBBase* a_gdb, int src_lev, int lev);

Vector<int> computeNeighborProcs (const ParGDBBase* a_gdb, int ngrow);

}

#endif // include guard
