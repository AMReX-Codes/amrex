#ifndef AMREX_PARTICLEUTIL_H_
#define AMREX_PARTICLEUTIL_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_MakeParticle.H>
#include <AMReX_Math.H>
#include <AMReX_MFIter.H>
#include <AMReX_ParGDB.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParticleBufferMap.H>
#include <AMReX_TypeTraits.H>
#include <AMReX_Scan.H>

#include <limits>

namespace amrex {

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
template <class Iterator, std::enable_if_t<IsParticleIterator<Iterator>::value && !Iterator::ContainerType::ParticleType::is_soa_particle, int> foo = 0>
int
numParticlesOutOfRange (Iterator const& pti, IntVect nGrow)
{
    using ParticleType = typename Iterator::ContainerType::ParticleType;

    const auto& tile = pti.GetParticleTile();
    const auto np = tile.numParticles();
    const auto& ptd = tile.getConstParticleTileData();
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
        auto p = make_particle<ParticleType>{}(ptd,i);
        if ((p.id() < 0)) { return false; }
        using AssignorType = typename Iterator::CellAssignor;
        AssignorType assignor;
        IntVect iv = assignor(p, plo, dxi, domain);
        return !box.contains(iv);
    });
    int hv = amrex::get<0>(reduce_data.value(reduce_op));
    return hv;
}

template <class Iterator, std::enable_if_t<IsParticleIterator<Iterator>::value && Iterator::ContainerType::ParticleType::is_soa_particle, int> foo = 0>
int
numParticlesOutOfRange (Iterator const& pti, IntVect nGrow)
{
    using ParticleType = typename Iterator::ContainerType::ConstParticleType;

    const auto& tile = pti.GetParticleTile();
    const auto tile_data = tile.getConstParticleTileData();
    const auto np = tile.numParticles();
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
        ParticleType p(tile_data,i);
        if ((p.id() < 0)) { return false; }
        using AssignorType = typename Iterator::CellAssignor;
        AssignorType assignor;
        IntVect iv = assignor(p, plo, dxi, domain);
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
    unsigned int operator() (const T& ptd, int i) const
    {
        auto p = ptd[i];
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

template <typename PTD>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
IntVect getParticleCell (PTD const& ptd, int i,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                         amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                         const Box& domain) noexcept
{
    if constexpr (PTD::ParticleType::is_soa_particle)
    {
        IntVect iv(
                   AMREX_D_DECL(int(amrex::Math::floor((ptd.m_rdata[0][i]-plo[0])*dxi[0])),
                                int(amrex::Math::floor((ptd.m_rdata[1][i]-plo[1])*dxi[1])),
                                int(amrex::Math::floor((ptd.m_rdata[2][i]-plo[2])*dxi[2]))));
        iv += domain.smallEnd();
        return iv;
    } else {
        return getParticleCell(ptd.m_aos[i], plo, dxi, domain);;
    }
}

struct DefaultAssignor
{

    template <typename P>
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (P const& p,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                        const Box& domain) const noexcept
    {
        return getParticleCell(p, plo, dxi, domain);
    }
};

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
int getParticleGrid (P const& p, amrex::Array4<int> const& mask,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                     amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                     const Box& domain) noexcept
{
    if (p.id() < 0) { return -1; }
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
        if (! is_per[idim]) { continue; }
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

/**
 * \brief Reorders the ParticleTile into two partitions
 * left [0, num_left-1] and right [num_left, ptile.numParticles()-1]
 * and returns the number of particles in the left partition.
 *
 * The functor is_left [(ParticleTileData ptd, int index) -> bool] maps each particle to
 * either the left [return true] or the right [return false] partition.
 * It must return the same result if evaluated multiple times for the same particle.
 *
 * \param ptile the ParticleTile to partition
 * \param is_left functor to map particles to a partition
 */
template <typename PTile, typename ParFunc>
int
partitionParticles (PTile& ptile, ParFunc const& is_left)
{
    const int np = ptile.numParticles();
    if (np == 0) { return 0; }

    auto ptd = ptile.getParticleTileData();

    const int num_left = Reduce::Sum<int>(np,
        [=] AMREX_GPU_DEVICE (int i) -> int
        {
            return int(is_left(ptd, i));
        });

    // The ptile will be partitioned into left [0, num_left-1] and right [num_left, np-1].
    //
    // Note that currently the number of particles in [0, num_left-1] that should belong to the
    // right partition is equal to the number of particles in [num_left, np-1] that should belong
    // in the left partition. We will define num_swaps to be this number. This is the minimum
    // number of swaps that need to be performed to partition the ptile in place for any algorithm.
    //
    // From this it is easy to see that
    // max_num_swaps = min(size([0, num_left-1]), size([num_left, np-1]))
    // is an upper bound for num_swaps.

    const int max_num_swaps = std::min(num_left, np - num_left);
    if (max_num_swaps == 0) { return num_left; }

    Gpu::DeviceVector<int> index_left(max_num_swaps);
    Gpu::DeviceVector<int> index_right(max_num_swaps);
    int * const p_index_left = index_left.dataPtr();
    int * const p_index_right = index_right.dataPtr();

    // The num_swaps particles that are in [0, num_left-1] but should be moved to the right
    // partition are at the same time the first num_swaps particles for which is_left is false
    // in all the ptile.
    // Similarly, the num_swaps particles in [num_left, np-1] that should be moved to the left
    // partition are the last num_swaps particles of the ptile for which is_left is true.
    //
    // The PrefixSum is used to find exactly these particles and store their indexes in
    // index_left and index_right. Since num_swaps is not known, the first max_num_swaps
    // particles are stored instead. Here, dst = num_left-1-(i-s) is used to effectively reverse
    // the PrefixSum to store the last particles for which is_left is true.
    //
    // This way, all indexes in index_right are in ascending order, and all indexes in
    // index_left are in descending order.

    Scan::PrefixSum<int>(np,
        [=] AMREX_GPU_DEVICE (int i) -> int
        {
            return int(!is_left(ptd, i));
        },
        [=] AMREX_GPU_DEVICE (int i, int const& s)
        {
            if (!is_left(ptd, i)) {
                int dst = s;
                if (dst < max_num_swaps) {
                    p_index_right[dst] = i;
                }
            } else {
                int dst = num_left-1-(i-s); // avoid integer overflow
                if (dst < max_num_swaps) {
                    p_index_left[dst] = i;
                }
            }
        },
        Scan::Type::exclusive, Scan::noRetSum);

    // Finally, the particles are swapped. Since max_num_swaps is only an upper bound for num_swaps,
    // some swaps should not be performed and need to be skipped. This is the case if the index
    // in index_left[i] is already in the left partition or the index in index_right[i] is already
    // in the right partition. These two cases coincide for the same i because index_right is in
    // ascending order and index_left in descending order. This means for both index_left and
    // index_right the first num_swaps particles need to be swapped, and the particles after that
    // should be skipped.
    //
    // The check right_i < left_i makes sure that the particle going to the right partition is
    // actually coming from the left partition, which has a lower index than the other particle and
    // visa-versa.
    //
    // Since exactly num_swaps swap operations are performed in the end, which is the smallest
    // number possible, this algorithm is optimal in the number of swap operations.
    // This results in good performance in practice if the size of a particle is large enough that
    // it compensates for the extra kernel launches and evaluations of is_left which this
    // algorithm needs.

    ParallelFor(max_num_swaps,
        [=] AMREX_GPU_DEVICE (int i)
        {
            int left_i = p_index_left[i];
            int right_i = p_index_right[i];
            if (right_i < left_i) {
                swapParticle(ptd, ptd, left_i, right_i);
            }
        });

    Gpu::streamSynchronize(); // for index_left and index_right deallocation

    return num_left;
}

template <typename PTile>
void
removeInvalidParticles (PTile& ptile)
{
    const int new_size = partitionParticles(ptile,
        [=] AMREX_GPU_DEVICE (auto& ptd, int i) {
            return ptd.id(i).is_valid();
        });
    ptile.resize(new_size);
}

#if defined(AMREX_USE_GPU)

template <typename PTile, typename PLocator, typename CellAssignor>
int
partitionParticlesByDest (PTile& ptile, const PLocator& ploc, CellAssignor const& assignor,
                          const ParticleBufferMap& pmap,
                          const GpuArray<Real,AMREX_SPACEDIM>& plo,
                          const GpuArray<Real,AMREX_SPACEDIM>& phi,
                          const GpuArray<ParticleReal,AMREX_SPACEDIM>& rlo,
                          const GpuArray<ParticleReal,AMREX_SPACEDIM>& rhi,
                          const GpuArray<int ,AMREX_SPACEDIM>& is_per,
                          int lev, int gid, int /*tid*/,
                          int lev_min, int lev_max, int nGrow, bool remove_negative)
{
    auto getPID = pmap.getPIDFunctor();
    int pid = ParallelContext::MyProcSub();

    Gpu::DeviceVector<uint8_t> particle_stays(ptile.numParticles());
    uint8_t * const p_particle_stays = particle_stays.dataPtr();
    auto ptd = ptile.getParticleTileData();

    // the function for determining if a particle stays on this grid is very slow,
    // so we cache it in particle_stays to avoid evaluating it multiple times.
    ParallelFor(ptile.numParticles(),
        [=] AMREX_GPU_DEVICE (int i)
        {
            int assigned_grid;
            int assigned_lev;

            if (!ptd.id(i).is_valid())
            {
                assigned_grid = -1;
                assigned_lev  = -1;
            }
            else
            {
                auto p_prime = ptd.getSuperParticle(i);
                enforcePeriodic(p_prime, plo, phi, rlo, rhi, is_per);
                auto tup_prime = ploc(p_prime, lev_min, lev_max, nGrow, assignor);
                assigned_grid = amrex::get<0>(tup_prime);
                assigned_lev  = amrex::get<1>(tup_prime);
                if (assigned_grid >= 0)
                {
                    AMREX_D_TERM(ptd.pos(0, i) = p_prime.pos(0);,
                                 ptd.pos(1, i) = p_prime.pos(1);,
                                 ptd.pos(2, i) = p_prime.pos(2););
                }
                else if (lev_min > 0)
                {
                    AMREX_D_TERM(p_prime.pos(0) = ptd.pos(0, i);,
                                 p_prime.pos(1) = ptd.pos(1, i);,
                                 p_prime.pos(2) = ptd.pos(2, i););
                    auto tup = ploc(p_prime, lev_min, lev_max, nGrow, assignor);
                    assigned_grid = amrex::get<0>(tup);
                    assigned_lev  = amrex::get<1>(tup);
                }
            }

            p_particle_stays[i] = uint8_t(
                ((remove_negative == false) && (!ptd.id(i).is_valid())) ||
                ((assigned_grid == gid) && (assigned_lev == lev) && (getPID(lev, gid) == pid)));
        });

    return partitionParticles(ptile,
        [=] AMREX_GPU_DEVICE (auto& /* ptd */, int i) -> bool {
            return p_particle_stays[i];
        });
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
                               index_type nbins, F const& f)
{
    BL_PROFILE("PermutationForDeposition()");

#if defined(AMREX_USE_HIP)
    // MI250X has a small L2 cache and is more tolerant of atomic add contention,
    // so we use a small block size of 64 and the compressed layout.
    static constexpr index_type gpu_block_size = 64;
    static constexpr bool compressed_layout = true;
#else
    // A100 has a larger L2 cache and is very sensitive to atomic add contention,
    // so we use a large bock size of 1024 and not the compressed layout.
    static constexpr index_type gpu_block_size = 1024;
    static constexpr bool compressed_layout = false;
#endif

    static constexpr index_type gpu_block_size_m1 = gpu_block_size - 1;
    static constexpr index_type llist_guard = std::numeric_limits<index_type>::max();

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
        [pllist_start,pllist_next,pperm,pglobal_idx] AMREX_GPU_DEVICE () {
            __shared__ index_type sdata[gpu_block_size];
            __shared__ index_type global_idx_start;
            __shared__ index_type idx_start;

            index_type current_idx = 0;

            if constexpr (compressed_layout) {
                // Compressed layout: subsequent sweeps of up to gpu_block_size contiguous particles
                // are put right next to each other, while without the compressed layout,
                // there can be other particle sweeps from different locations between them.
                current_idx = pllist_start[threadIdx.x + gpu_block_size * blockIdx.x];

                index_type num_particles_thread = 0;
                while (current_idx != llist_guard) {
                    ++num_particles_thread;
                    current_idx = pllist_next[current_idx];
                }

                index_type num_particles_block =
                    Gpu::blockReduceSum<gpu_block_size>(num_particles_thread);

                if (threadIdx.x == 0) {
                    global_idx_start = Gpu::Atomic::Add(pglobal_idx, num_particles_block);
                }
            }

            current_idx = pllist_start[threadIdx.x + gpu_block_size * blockIdx.x];

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
                if (threadIdx.x == gpu_block_size_m1) {
                    if constexpr (compressed_layout) {
                        idx_start = global_idx_start;
                        global_idx_start += sdata[gpu_block_size_m1];
                    } else {
                        idx_start = Gpu::Atomic::Add(pglobal_idx, sdata[gpu_block_size_m1]);
                    }
                }
                __syncthreads();
                sdata[threadIdx.x] += idx_start;
                if (current_idx != llist_guard) {
                    pperm[sdata[threadIdx.x] - 1] = current_idx;
                    current_idx = pllist_next[current_idx];
                }
            }
        });
#else
    amrex::ignore_unused(pperm, pglobal_idx, compressed_layout);
    Abort("PermutationForDeposition only implemented for CUDA and HIP");
#endif

    Gpu::Device::streamSynchronize();
}

template <class index_type, class PTile>
void PermutationForDeposition (Gpu::DeviceVector<index_type>& perm, index_type nitems,
                               const PTile& ptile, Box bx, Geometry geom, const IntVect idx_type)
{
    AMREX_ALWAYS_ASSERT(idx_type.allGE(0) && idx_type.allLE(2));

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

    auto ptd = ptile.getConstParticleTileData();
    using ParticleType = typename PTile::ParticleType::ConstType;
    PermutationForDeposition<index_type>(perm, nitems, bx.numPts() * ref_product,
        [=] AMREX_GPU_DEVICE (index_type idx) noexcept
            {
                const auto& p = make_particle<ParticleType>{}(ptd,idx);

                IntVect iv = ((p.pos() - pos_offset) * dxi).round();

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
