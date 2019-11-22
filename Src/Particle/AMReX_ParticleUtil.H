#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_TypeTraits.H>

#include <limits>

namespace amrex
{

/**
 * \brief Returns the number of particles that are more than nGrow cells
 * from the box correspond to the input iterator.
 *
 * \tparam Iterator an AMReX ParticleIterator
 *
 * \param the iterator pointing to the current grid/tile to test
 * \param nGrow the number of grow cells allowed.
 *
 */
template <class Iterator, EnableIf_t<amrex::IsParticleIterator<Iterator>::value, int> foo = 0>
int
numParticlesOutOfRange (Iterator const& pti, int nGrow)
{
    using ParticleType = typename Iterator::ContainerType::ParticleType;

    const auto& tile = pti.GetParticleTile();
    const auto np = tile.numParticles();
    const auto& aos = tile.GetArrayOfStructs();
    const auto pstruct = aos().dataPtr();
    const auto& geom = pti.Geom(pti.GetLevel());

    const auto domain = geom.Domain();
    const auto plo = geom.ProbLoArray();
    const auto dxi = geom.InvCellSizeArray();
    
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
    return hv;
}
    
/**
 * \brief Returns the number of particles that are more than nGrow cells
 * from their assigned box.
 *
 * This version tests over all levels.
 *
 * \tparam PC a type of AMReX particle container.
 *
 * \param pc the particle container to test
 * \param nGrow the number of grow cells allowed.
 *
 */    
template <class PC, EnableIf_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
int
numParticlesOutOfRange (PC const& pc, int nGrow)
{
    return numParticlesOutOfRange(pc, 0, pc.finestLevel(), nGrow);
}

/**
 * \brief Returns the number of particles that are more than nGrow cells
 * from their assigned box.
 *
 * This version goes over only the specified levels
 *
 * \tparam PC a type of AMReX particle container.
 *
 * \param pc the particle container to test
 * \param lev_min the minimum level to test
 * \param lev_max the maximum level to test
 * \param nGrow the number of grow cells allowed.
 *
 */        
template <class PC, EnableIf_t<amrex::IsParticleContainer<PC>::value, int> foo = 0>
int
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
            num_wrong += numParticlesOutOfRange(pti, nGrow);
        }
    }
    ParallelDescriptor::ReduceIntSum(num_wrong);

    return num_wrong;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getTileIndex (const IntVect& iv, const Box& box, const bool a_do_tiling, 
		  const IntVect& a_tile_size, Box& tbx)
{
    if (a_do_tiling == false) {
        tbx = box;
        return 0;
    } else {
        //
        // This function must be consistent with FabArrayBase::buildTileArray function!!!
        //
        auto tiling_1d = [](int i, int lo, int hi, int tilesize,
                            int& ntile, int& tileidx, int& tlo, int& thi) {
            int ncells = hi-lo+1;
            ntile = amrex::max(ncells/tilesize, 1);
            int ts_right = ncells/ntile;
            int ts_left  = ts_right+1;
            int nleft = ncells - ntile*ts_right;
	    int ii = i - lo;
            int nbndry = nleft*ts_left;
            if (ii < nbndry) {
                tileidx = ii / ts_left; // tiles on the left of nbndry have size of ts_left
                tlo = lo + tileidx * ts_left;
                thi = tlo + ts_left - 1;
            } else {
                tileidx = nleft + (ii-nbndry) / ts_right;  // tiles on the right: ts_right
                tlo = lo + tileidx * ts_right + nleft;
                thi = tlo + ts_right - 1;
            }
        };
        const IntVect& small = box.smallEnd();
        const IntVect& big   = box.bigEnd();
        IntVect ntiles, ivIndex, tilelo, tilehi;

        AMREX_D_TERM(int iv0 = amrex::min(amrex::max(iv[0], small[0]), big[0]);,
		     int iv1 = amrex::min(amrex::max(iv[1], small[1]), big[1]);,
		     int iv2 = amrex::min(amrex::max(iv[2], small[2]), big[2]););

        AMREX_D_TERM(tiling_1d(iv0, small[0], big[0], a_tile_size[0], ntiles[0], ivIndex[0], tilelo[0], tilehi[0]);,
		     tiling_1d(iv1, small[1], big[1], a_tile_size[1], ntiles[1], ivIndex[1], tilelo[1], tilehi[1]);,
		     tiling_1d(iv2, small[2], big[2], a_tile_size[2], ntiles[2], ivIndex[2], tilelo[2], tilehi[2]););

        tbx = Box(tilelo, tilehi);

        return AMREX_D_TERM(ivIndex[0], + ntiles[0]*ivIndex[1], + ntiles[0]*ntiles[1]*ivIndex[2]);
    }
}

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

IntVect computeRefFac (const ParGDBBase* a_gdb, int src_lev, int lev);

Vector<int> computeNeighborProcs (const ParGDBBase* a_gdb, int ngrow);

}

#endif // include guard
