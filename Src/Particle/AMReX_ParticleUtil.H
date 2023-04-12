#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Scan.H>

#include <limits>

namespace amrex
{

namespace particle_detail {

template <typename F, typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f, P const& p, amrex::RandomEngine const& engine) noexcept
    -> decltype(f(P{},RandomEngine{}))
{
    return f(p,engine);
}

template <typename F, typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f, P const& p, amrex::RandomEngine const&) noexcept
    -> decltype(f(P{}))
{
    return f(p);
}

template <typename F, typename SrcData, typename N>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f, SrcData const& src, N i, amrex::RandomEngine const& engine) noexcept
    -> decltype(f(SrcData{},N{},RandomEngine{}))
{
    return f(src,i,engine);
}

template <typename F, typename SrcData, typename N>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f, SrcData const& src, N i, amrex::RandomEngine const&) noexcept
    -> decltype(f(SrcData{},N{}))
{
    return f(src,i);
}

// The next several functions are used by ParticleReduce

// Lambda takes a Particle
template <typename F, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i) noexcept
    -> decltype(f(p.m_aos[i]))
{
    return f(p.m_aos[i]);
}

// Lambda takes a SuperParticle
template <typename F, int NSR, int NSI, int NAR, int NAI,
          typename std::enable_if<NAR != 0 || NAI != 0, int>::type = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i) noexcept
    -> decltype(f(p.getSuperParticle(i)))
{
    return f(p.getSuperParticle(i));
}

// Lambda takes a ConstParticleTileData
template <typename F, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i) noexcept
    -> decltype(f(p, i))
{
    return f(p, i);
}

// These next several functions are used by ParticleToMesh and MeshToParticle

// Lambda takes a Particle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
    -> decltype(f(p.m_aos[i], fabarr, plo, dxi))
{
    return f(p.m_aos[i], fabarr, plo, dxi);
}

// Lambda takes a Particle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const&,
             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
    -> decltype(f(p.m_aos[i], fabarr))
{
    return f(p.m_aos[i], fabarr);
}

// Lambda takes a Particle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<const T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
    -> decltype(f(p.m_aos[i], fabarr, plo, dxi))
{
    return f(p.m_aos[i], fabarr, plo, dxi);
}

// Lambda takes a Particle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<const T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const&,
             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
    -> decltype(f(p.m_aos[i], fabarr))
{
    return f(p.m_aos[i], fabarr);
}

// Lambda takes a SuperParticle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI,
          typename std::enable_if<(NAR != 0) || (NAI != 0), int>::type = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
    -> decltype(f(p.getSuperParticle(i), fabarr, plo, dxi))
{
    return f(p.getSuperParticle(i), fabarr, plo, dxi);
}

// Lambda takes a SuperParticle
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI,
          typename std::enable_if<(NAR != 0) || (NAI != 0), int>::type = 0>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const&,
             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
    -> decltype(f(p.getSuperParticle(i), fabarr))
{
    return f(p.getSuperParticle(i), fabarr);
}

// Lambda takes a ConstParticleTileData
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const&,
             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
    -> decltype(f(p, i, fabarr))
{
    return f(p, i, fabarr);
}

// Lambda takes a ConstParticleTileData
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ConstParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
    -> decltype(f(p, i, fabarr, plo, dxi))
{
    return f(p, i, fabarr, plo, dxi);
}

// Lambda takes a ParticleTileData
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<const T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const& plo,
             GpuArray<Real,AMREX_SPACEDIM> const& dxi) noexcept
    -> decltype(f(p, i, fabarr, plo, dxi))
{
    return f(p, i, fabarr, plo, dxi);
}

// Lambda takes a ParticleTileData
template <typename F, typename T, int NSR, int NSI, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
auto call_f (F const& f,
             const ParticleTileData<NSR, NSI, NAR, NAI>& p,
             const int i, Array4<const T> const& fabarr,
             GpuArray<Real,AMREX_SPACEDIM> const&,
             GpuArray<Real,AMREX_SPACEDIM> const&) noexcept
    -> decltype(f(p, i, fabarr))
{
    return f(p, i, fabarr);
}


}

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
template <class Iterator, std::enable_if_t<IsParticleIterator<Iterator>::value, int> foo = 0>
int
numParticlesOutOfRange (Iterator const& pti, int nGrow)
{
    return numParticlesOutOfRange(pti,
                                  IntVect(AMREX_D_DECL(nGrow, nGrow, nGrow)));
}

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
template <class Iterator, std::enable_if_t<IsParticleIterator<Iterator>::value, int> foo = 0>
int
numParticlesOutOfRange (Iterator const& pti, IntVect nGrow)
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
        IntVect iv(
            AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                         int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                         int(amrex::Math::floor((p.pos(2)-plo[2])*dxi[2]))));
        iv += domain.smallEnd();
        return !box.contains(iv);
    });
    int hv = amrex::get<0>(reduce_data.value(reduce_op));
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
template <class PC, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
int
numParticlesOutOfRange (PC const& pc, int nGrow)
{
    return numParticlesOutOfRange(pc, 0, pc.finestLevel(), nGrow);
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
template <class PC, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
int
numParticlesOutOfRange (PC const& pc, IntVect nGrow)
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
template <class PC, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
int
numParticlesOutOfRange (PC const& pc, int lev_min, int lev_max, int nGrow)
{
    BL_PROFILE("numParticlesOutOfRange()");

    return numParticlesOutOfRange(pc, lev_min, lev_max,
                                  IntVect(AMREX_D_DECL(nGrow, nGrow, nGrow)));
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
template <class PC, std::enable_if_t<IsParticleContainer<PC>::value, int> foo = 0>
int
numParticlesOutOfRange (PC const& pc, int lev_min, int lev_max, IntVect nGrow)
{
    BL_PROFILE("numParticlesOutOfRange()");

    using ParIter = typename PC::ParConstIterType;
    int num_wrong = 0;
    for (int lev = lev_min; lev <= lev_max; ++lev)
    {
#ifdef AMREX_USE_OMP
#pragma omp parallel if (Gpu::notInLaunchRegion() && !system::regtest_reduction) reduction(+:num_wrong)
#endif
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            num_wrong += numParticlesOutOfRange(pti, nGrow);
        }
    }
    ParallelAllReduce::Sum(num_wrong, ParallelContext::CommunicatorSub());

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

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int numTilesInBox (const Box& box, const bool a_do_tiling, const IntVect& a_tile_size)
{
    if (a_do_tiling == false) {
        return 1;
    } else {
        //
        // This function must be consistent with FabArrayBase::buildTileArray function!!!
        //
        auto tiling_1d = [](int lo, int hi, int tilesize, int& ntile) {
            int ncells = hi-lo+1;
            ntile = amrex::max(ncells/tilesize, 1);
        };

        const IntVect& small = box.smallEnd();
        const IntVect& big   = box.bigEnd();
        IntVect ntiles;

        AMREX_D_TERM(tiling_1d(small[0], big[0], a_tile_size[0], ntiles[0]);,
                     tiling_1d(small[1], big[1], a_tile_size[1], ntiles[1]);,
                     tiling_1d(small[2], big[2], a_tile_size[2], ntiles[2]););

        return AMREX_D_TERM(ntiles[0], *=ntiles[1], *=ntiles[2]);
    }
}

struct BinMapper
{
    BinMapper(const int* off_bins_p,
              const GpuArray<Real,AMREX_SPACEDIM>* dxi_p,
              const GpuArray<Real,AMREX_SPACEDIM>* plo_p,
              const Dim3* lo_p,
              const Dim3* hi_p,
              int* bin_type_array=nullptr)
        : m_off_bins_p(off_bins_p), m_dxi_p(dxi_p), m_plo_p(plo_p)                  ,
          m_lo_p(lo_p)            , m_hi_p(hi_p)  , m_bin_type_array(bin_type_array) {}

    template <typename T>
    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (const T& p, int i) const
    {
        int type   = (m_bin_type_array) ? m_bin_type_array[i] : 0;
        int offset = m_off_bins_p[type];

        AMREX_D_TERM(AMREX_ASSERT((p.pos(0)-m_plo_p[type][0])*m_dxi_p[type][0] - m_lo_p[type].x >= 0.0);,
                     AMREX_ASSERT((p.pos(1)-m_plo_p[type][1])*m_dxi_p[type][1] - m_lo_p[type].y >= 0.0);,
                     AMREX_ASSERT((p.pos(2)-m_plo_p[type][2])*m_dxi_p[type][2] - m_lo_p[type].z >= 0.0));

        auto iv = IntVect(AMREX_D_DECL(static_cast<int>(amrex::Math::floor((p.pos(0)-m_plo_p[type][0])*m_dxi_p[type][0])) - m_lo_p[type].x,
                                       static_cast<int>(amrex::Math::floor((p.pos(1)-m_plo_p[type][1])*m_dxi_p[type][1])) - m_lo_p[type].y,
                                       static_cast<int>(amrex::Math::floor((p.pos(2)-m_plo_p[type][2])*m_dxi_p[type][2])) - m_lo_p[type].z));
        auto iv3 = iv.dim3();
        int nx   = m_hi_p[type].x-m_lo_p[type].x+1;
        int ny   = m_hi_p[type].y-m_lo_p[type].y+1;
        int nz   = m_hi_p[type].z-m_lo_p[type].z+1;
        int uix = amrex::min(nx-1,amrex::max(0,iv3.x));
        int uiy = amrex::min(ny-1,amrex::max(0,iv3.y));
        int uiz = amrex::min(nz-1,amrex::max(0,iv3.z));
        return static_cast<unsigned int>( (uix * ny + uiy) * nz + uiz + offset );
    }

private:
    const int* m_off_bins_p;
    const GpuArray<Real,AMREX_SPACEDIM>* m_dxi_p;
    const GpuArray<Real,AMREX_SPACEDIM>* m_plo_p;
    const Dim3* m_lo_p;
    const Dim3* m_hi_p;
    int* m_bin_type_array;
};

struct GetParticleBin
{
    GpuArray<Real,AMREX_SPACEDIM> plo;
    GpuArray<Real,AMREX_SPACEDIM> dxi;
    Box domain;
    IntVect bin_size;
    Box box;

    template <typename ParticleType>
    AMREX_GPU_HOST_DEVICE
    unsigned int operator() (const ParticleType& p) const noexcept
    {
        Box tbx;
        auto iv = getParticleCell(p, plo, dxi, domain);
        auto tid = getTileIndex(iv, box, true, bin_size, tbx);
        return static_cast<unsigned int>(tid);
    }
};

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getParticleCell (P const& p,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    IntVect iv(
        AMREX_D_DECL(int(amrex::Math::floor((p.pos(0)-plo[0])*dxi[0])),
                     int(amrex::Math::floor((p.pos(1)-plo[1])*dxi[1])),
                     int(amrex::Math::floor((p.pos(2)-plo[2])*dxi[2]))));
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
bool enforcePeriodic (P& p,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                      amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& phi,
                      amrex::GpuArray<amrex::ParticleReal,AMREX_SPACEDIM> const& rlo,
                      amrex::GpuArray<amrex::ParticleReal,AMREX_SPACEDIM> const& rhi,
                      amrex::GpuArray<int,AMREX_SPACEDIM> const& is_per) noexcept
{
    bool shifted = false;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim)
    {
        if (! is_per[idim]) continue;
        if (p.pos(idim) > rhi[idim]) {
            while (p.pos(idim) > rhi[idim]) {
                p.pos(idim) -= static_cast<ParticleReal>(phi[idim] - plo[idim]);
            }
            // clamp to avoid precision issues;
            if (p.pos(idim) < rlo[idim]) {
                p.pos(idim) = rlo[idim];
            }
            shifted = true;
        }
        else if (p.pos(idim) < rlo[idim]) {
            while (p.pos(idim) < rlo[idim]) {
                p.pos(idim) += static_cast<ParticleReal>(phi[idim] - plo[idim]);
            }
            // clamp to avoid precision issues;
            if (p.pos(idim) > rhi[idim]) {
                p.pos(idim) = rhi[idim];
            }
            shifted = true;
        }
        AMREX_ASSERT( (p.pos(idim) >= rlo[idim] ) && ( p.pos(idim) <= rhi[idim] ));
    }

    return shifted;
}

#if defined(AMREX_USE_GPU)

template <typename PTile, typename PLocator>
int
partitionParticlesByDest (PTile& ptile, const PLocator& ploc, const ParticleBufferMap& pmap,
                          const GpuArray<Real,AMREX_SPACEDIM>& plo,
                          const GpuArray<Real,AMREX_SPACEDIM>& phi,
                          const GpuArray<ParticleReal,AMREX_SPACEDIM>& rlo,
                          const GpuArray<ParticleReal,AMREX_SPACEDIM>& rhi,
                          const GpuArray<int ,AMREX_SPACEDIM>& is_per,
                          int lev, int gid, int /*tid*/,
                          int lev_min, int lev_max, int nGrow, bool remove_negative)
{
    auto& aos = ptile.GetArrayOfStructs();
    const int np = aos.numParticles();

    if (np == 0) return 0;

    auto getPID = pmap.getPIDFunctor();
    auto p_ptr = &(aos[0]);

    int pid = ParallelContext::MyProcSub();
    constexpr int chunk_size = 256*256*256;
    int num_chunks = std::max(1, (np + (chunk_size - 1)) / chunk_size);

    PTile ptile_tmp;
    ptile_tmp.define(ptile.NumRuntimeRealComps(), ptile.NumRuntimeIntComps());
    ptile_tmp.resize(std::min(np, chunk_size));

    auto src_data = ptile.getParticleTileData();
    auto dst_data = ptile_tmp.getParticleTileData();

    int last_offset = 0;
    for (int ichunk = 0; ichunk < num_chunks; ++ichunk)
    {
        int this_offset = ichunk*chunk_size;
        int this_chunk_size = std::min(chunk_size, np - this_offset);

        int num_stay;
        {
            auto particle_stays = [=] AMREX_GPU_DEVICE (int i) -> int
            {
                int assigned_grid;
                int assigned_lev;

                auto& p = p_ptr[i+this_offset];

                if (p.id() < 0 )
                {
                    assigned_grid = -1;
                    assigned_lev  = -1;
                }
                else
                {
                    auto p_prime = p;
                    enforcePeriodic(p_prime, plo, phi, rlo, rhi, is_per);
                    auto tup_prime = ploc(p_prime, lev_min, lev_max, nGrow);
                    assigned_grid = amrex::get<0>(tup_prime);
                    assigned_lev  = amrex::get<1>(tup_prime);
                    if (assigned_grid >= 0)
                    {
                      AMREX_D_TERM(p.pos(0) = p_prime.pos(0);,
                                   p.pos(1) = p_prime.pos(1);,
                                   p.pos(2) = p_prime.pos(2););
                    }
                    else if (lev_min > 0)
                    {
                      auto tup = ploc(p, lev_min, lev_max, nGrow);
                      assigned_grid = amrex::get<0>(tup);
                      assigned_lev  = amrex::get<1>(tup);
                    }
                }

                if ((remove_negative == false) && (p.id() < 0)) {
                    return true;
                }

                return ((assigned_grid == gid) && (assigned_lev == lev) && (getPID(lev, gid) == pid));
            };

            num_stay = Scan::PrefixSum<int> (this_chunk_size,
                          [=] AMREX_GPU_DEVICE (int i) -> int
                          {
                              return particle_stays(i);
                          },
                          [=] AMREX_GPU_DEVICE (int i, int const& s)
                          {
                              int src_i = i + this_offset;
                              int dst_i = particle_stays(i) ? s : this_chunk_size-1-(i-s);
                              copyParticle(dst_data, src_data, src_i, dst_i);
                          },
                          Scan::Type::exclusive);
        }

        if (num_chunks == 1)
        {
            ptile.swap(ptile_tmp);
        }
        else
        {
            AMREX_FOR_1D(this_chunk_size, i,
                         {
                             copyParticle(src_data, dst_data, i, i + this_offset);
                         });
        }

        if ( ichunk > 0 )
        {
            int num_swap = std::min(this_offset - last_offset, num_stay);
            AMREX_FOR_1D( num_swap, i,
            {
                swapParticle(src_data, src_data, last_offset + i,
                             this_offset + num_stay - 1 - i);
            });
        }

        last_offset += num_stay;
    }

    return last_offset;
}

#endif

template <class PC1, class PC2>
bool SameIteratorsOK (const PC1& pc1, const PC2& pc2) {
    if (pc1.numLevels() != pc2.numLevels()) {return false;}
    if (pc1.do_tiling != pc2.do_tiling) {return false;}
    if (pc1.tile_size != pc2.tile_size) {return false;}
    for (int lev = 0; lev < pc1.numLevels(); ++lev) {
        if (pc1.ParticleBoxArray(lev) != pc2.ParticleBoxArray(lev)) {return false;}
        if (pc1.ParticleDistributionMap(lev) != pc2.ParticleDistributionMap(lev)) {return false;}
    }
    return true;
}

template <class PC>
void EnsureThreadSafeTiles(PC& pc) {
    using Iter = typename PC::ParIterType;
    for (int lev = 0; lev < pc.numLevels(); ++lev) {
        for (Iter pti(pc, lev); pti.isValid(); ++pti) {
            pc.DefineAndReturnParticleTile(lev, pti);
        }
    }
}

IntVect computeRefFac (const ParGDBBase* a_gdb, int src_lev, int lev);

Vector<int> computeNeighborProcs (const ParGDBBase* a_gdb, int ngrow);

namespace particle_detail
{
template <typename C>
void clearEmptyEntries (C& c)
{
    for (auto c_it = c.begin(); c_it != c.end(); /* no ++ */)
    {
        if (c_it->second.empty()) { c.erase(c_it++); }
        else { ++c_it; }
    }
}
}

template <class index_type, typename F>
void PermutationForDeposition (Gpu::DeviceVector<index_type>& perm, index_type nitems,
                               index_type nbins, F&& f)
{
    BL_PROFILE("PermutationForDeposition()");

    constexpr index_type gpu_block_size = 1024;
    constexpr index_type gpu_block_size_m1 = gpu_block_size - 1;
    constexpr index_type llist_guard = std::numeric_limits<index_type>::max();

    // round up to gpu_block_size
    nbins = (nbins + gpu_block_size_m1) / gpu_block_size * gpu_block_size;

    Gpu::DeviceVector<index_type> llist_start(nbins, llist_guard);
    Gpu::DeviceVector<index_type> llist_next(nitems);
    perm.resize(nitems);
    Gpu::DeviceScalar<index_type> global_idx(0);

    index_type* pllist_start = llist_start.dataPtr();
    index_type* pllist_next = llist_next.dataPtr();
    index_type* pperm = perm.dataPtr();
    index_type* pglobal_idx = global_idx.dataPtr();

    amrex::ParallelFor(nitems, [=] AMREX_GPU_DEVICE (index_type i) noexcept
    {
        i = nitems - i - 1;
        pllist_next[i] = Gpu::Atomic::Exch(pllist_start + f(i), i);
    });

#if defined(AMREX_USE_CUDA) || defined(AMREX_USE_HIP)
    amrex::launch<gpu_block_size>(nbins / gpu_block_size, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () {
            __shared__ index_type sdata[gpu_block_size];
            index_type current_idx = pllist_start[threadIdx.x + gpu_block_size * blockIdx.x];

            while (true) {
                sdata[threadIdx.x] = index_type(current_idx != llist_guard);
                index_type x = 0;

                // simple block wide prefix sum
                for (index_type i = 1; i<gpu_block_size; i*=2) {
                    __syncthreads();
                    if (threadIdx.x >= i) {
                        x = sdata[threadIdx.x - i];
                    }
                    __syncthreads();
                    if (threadIdx.x >= i) {
                        sdata[threadIdx.x] += x;
                    }
                }
                __syncthreads();
                if (sdata[gpu_block_size_m1] == 0) {
                    break;
                }
                __syncthreads();
                if (threadIdx.x == gpu_block_size_m1) {
                    x = sdata[gpu_block_size_m1];
                    sdata[gpu_block_size_m1] = Gpu::Atomic::Add(pglobal_idx, x);
                }
                __syncthreads();
                if (threadIdx.x < gpu_block_size_m1) {
                    sdata[threadIdx.x] += sdata[gpu_block_size_m1];
                }
                __syncthreads();
                if (threadIdx.x == gpu_block_size_m1) {
                    sdata[gpu_block_size_m1] += x;
                }
                __syncthreads();

                if (current_idx != llist_guard) {
                    pperm[sdata[threadIdx.x] - 1] = current_idx;
                    current_idx = pllist_next[current_idx];
                }
            }
        });
#else
    Abort("Not implemented");
#endif

    Gpu::Device::streamSynchronize();
}

template <class index_type, class pos_struct>
void PermutationForDeposition (Gpu::DeviceVector<index_type>& perm, index_type nitems,
                               pos_struct const* v, Box bx, Geometry geom, const IntVect idx_type)
{
    AMREX_ALWAYS_ASSERT(idx_type.allGE(IntVect(0)) && idx_type.allLE(IntVect(2)));

    const IntVect refine_vect = max(idx_type, IntVect(1)).min(IntVect(2));
    const IntVect type_vect = idx_type - idx_type / 2 * 2;

    geom.refine(refine_vect);

    Box domain = geom.Domain();

    bx.convert(type_vect);
    domain.convert(type_vect);

    const RealVect dxi(geom.InvCellSize());
    const RealVect pos_offset = Real(0.5) * (RealVect(geom.ProbLo()) + RealVect(geom.ProbHi())
        - RealVect(geom.CellSize()) * RealVect(domain.smallEnd() + domain.bigEnd()));

    const int ref_product = AMREX_D_TERM(refine_vect[0], * refine_vect[1], * refine_vect[2]);
    const IntVect ref_offset(AMREX_D_DECL(1, refine_vect[0], refine_vect[0] * refine_vect[1]));

    PermutationForDeposition<index_type>(perm, nitems, bx.numPts() * ref_product,
        [=] AMREX_GPU_DEVICE (index_type idx) noexcept
        {
            IntVect iv = ((v[idx].pos() - pos_offset) * dxi).round();

            IntVect iv_coarse = iv / refine_vect;
            IntVect iv_remainder = iv - iv_coarse * refine_vect;

            iv_coarse = iv_coarse.max(bx.smallEnd());
            iv_coarse = iv_coarse.min(bx.bigEnd());
            return bx.index(iv_coarse) + bx.numPts() * (iv_remainder * ref_offset).sum();
        });
}


#ifdef AMREX_USE_HDF5_ASYNC
void async_vol_es_wait_particle();
void async_vol_es_wait_close_particle();
#endif
}

#endif // include guard
