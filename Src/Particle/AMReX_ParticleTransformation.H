#ifndef AMREX_PARTICLETRANSFORMATION_H_
#define AMREX_PARTICLETRANSFORMATION_H_
#include <AMReX_Config.H>

#include <AMReX_IntVect.H>
#include <AMReX_Box.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParticleUtil.H>

namespace amrex
{

/**
 * \brief A general single particle copying routine that can run on the GPU.
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
template <typename T_ParticleType, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void copyParticle (const      ParticleTileData<T_ParticleType, NAR, NAI>& dst,
                   const ConstParticleTileData<T_ParticleType, NAR, NAI>& src,
                   int src_i, int dst_i) noexcept
{
    AMREX_ASSERT(dst.m_num_runtime_real == src.m_num_runtime_real);
    AMREX_ASSERT(dst.m_num_runtime_int  == src.m_num_runtime_int );

    if constexpr(!T_ParticleType::is_soa_particle) {
        dst.m_aos[dst_i] = src.m_aos[src_i];
    } else {
        dst.m_idcpu[dst_i] = src.m_idcpu[src_i];
    }
    if constexpr(NAR > 0) {
        for (int j = 0; j < NAR; ++j) {
            dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
        }
    }
    for (int j = 0; j < dst.m_num_runtime_real; ++j) {
        dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
    }
    if constexpr(NAI > 0) {
        for (int j = 0; j < NAI; ++j) {
            dst.m_idata[j][dst_i] = src.m_idata[j][src_i];
        }
    }
    for (int j = 0; j < dst.m_num_runtime_int; ++j) {
        dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
    }
}

/**
 * \brief A general single particle copying routine that can run on the GPU.
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
template <typename T_ParticleType, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void copyParticle (const ParticleTileData<T_ParticleType, NAR, NAI>& dst,
                   const ParticleTileData<T_ParticleType, NAR, NAI>& src,
                   int src_i, int dst_i) noexcept
{
    AMREX_ASSERT(dst.m_num_runtime_real == src.m_num_runtime_real);
    AMREX_ASSERT(dst.m_num_runtime_int  == src.m_num_runtime_int );

    if constexpr(T_ParticleType::is_soa_particle) {
        dst.m_idcpu[dst_i] = src.m_idcpu[src_i];
    } else {
        dst.m_aos[dst_i] = src.m_aos[src_i];
    }
    for (int j = 0; j < NAR; ++j) {
        dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
    }
    for (int j = 0; j < dst.m_num_runtime_real; ++j) {
        dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
    }
    for (int j = 0; j < NAI; ++j) {
        dst.m_idata[j][dst_i] = src.m_idata[j][src_i];
    }
    for (int j = 0; j < dst.m_num_runtime_int; ++j) {
        dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
    }
}

/**
 * \brief A general single particle swapping routine that can run on the GPU.
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
template <typename T_ParticleType, int NAR, int NAI>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void swapParticle (const ParticleTileData<T_ParticleType, NAR, NAI>& dst,
                   const ParticleTileData<T_ParticleType, NAR, NAI>& src,
                   int src_i, int dst_i) noexcept
{
    AMREX_ASSERT(dst.m_num_runtime_real == src.m_num_runtime_real);
    AMREX_ASSERT(dst.m_num_runtime_int  == src.m_num_runtime_int );

    if constexpr(T_ParticleType::is_soa_particle) {
        amrex::Swap(src.m_idcpu[src_i], dst.m_idcpu[dst_i]);
    } else {
        amrex::Swap(src.m_aos[src_i], dst.m_aos[dst_i]);
    }
    for (int j = 0; j < NAR; ++j) {
        amrex::Swap(dst.m_rdata[j][dst_i], src.m_rdata[j][src_i]);
    }
    for (int j = 0; j < dst.m_num_runtime_real; ++j) {
        amrex::Swap(dst.m_runtime_rdata[j][dst_i], src.m_runtime_rdata[j][src_i]);
    }
    for (int j = 0; j < NAI; ++j) {
        amrex::Swap(dst.m_idata[j][dst_i], src.m_idata[j][src_i]);
    }
    for (int j = 0; j < dst.m_num_runtime_int; ++j) {
        amrex::Swap(dst.m_runtime_idata[j][dst_i], src.m_runtime_idata[j][src_i]);
    }
}

/**
 * \brief Copy particles from src to dst. This version copies all the
 * particles, writing them to the beginning of dst.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */
template <typename DstTile, typename SrcTile>
void copyParticles (DstTile& dst, const SrcTile& src) noexcept
{
    auto np = src.numParticles();
    copyParticles(dst, src, 0, 0, np);
}

/**
 * \brief Copy particles from src to dst. This version copies n particles
 * starting at index src_start, writing the result starting at dst_start.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam N the size type, e.g. Long
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param n the number of particles to write
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename N,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
void copyParticles (DstTile& dst, const SrcTile& src,
                    Index src_start, Index dst_start, N n) noexcept
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        copyParticle(dst_data, src_data, src_start+i, dst_start+i);
    });

    Gpu::streamSynchronize();
}

/**
 * \brief Apply the function f to all the particles in src, writing the
 * result to dst. This version does all the particles in src.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */
template <typename DstTile, typename SrcTile, typename F>
void transformParticles (DstTile& dst, const SrcTile& src, F&& f) noexcept
{
    auto np = src.numParticles();
    transformParticles(dst, src, 0, 0, np, std::forward<F>(f));
}

/**
 * \brief Apply the function f to particles in src, writing the
 * result to dst. This version applies the function to n particles
 * starting at index src_start, writing the result starting at dst_start.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam N the size type, e.g. Long
 * \tparam F a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param f the function that will be applied to each particle
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename N, typename F,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
void transformParticles (DstTile& dst, const SrcTile& src,
                         Index src_start, Index dst_start, N n, F const& f) noexcept
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        f(dst_data, src_data, src_start+i, dst_start+i);
    });

    Gpu::streamSynchronize();
}

/**
 * \brief Apply the function f to all the particles in src, writing the
 * results to dst1 and dst2. This version does all the particles in src.
 *
 * \tparam DstTile1 the dst1 particle tile type
 * \tparam DstTile2 the dst2 particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam F a function object
 *
 * \param dst1 the first destination tile
 * \param dst2 the second destination tile
 * \param src the source tile
 * \param f the function that will be applied to each particle
 *
 */
template <typename DstTile1, typename DstTile2, typename SrcTile, typename F>
void transformParticles (DstTile1& dst1, DstTile2& dst2, const SrcTile& src, F&& f) noexcept
{
    auto np = src.numParticles();
    transformParticles(dst1, dst2, src, 0, 0, 0, np, std::forward<F>(f));
}

/**
 * \brief Apply the function f to particles in src, writing the
 * results to dst1 and dst2. This version applies the function to n particles
 * starting at index src_start, writing the result starting at dst1_start and dst2_start.
 *
 * \tparam DstTile1 the dst1 particle tile type
 * \tparam DstTile2 the dst2 particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam N the size type, e.g. Long
 * \tparam F a function object
 *
 * \param dst1 the first destination tile
 * \param dst2 the second destination tile
 * \param src the source tile
 * \param src_start the offset at which to start reading particles from src
 * \param dst1_start the offset at which to start writing particles to dst1
 * \param dst2_start the offset at which to start writing particles to dst2
 * \param f the function that will be applied to each particle
 *
 */
template <typename DstTile1, typename DstTile2, typename SrcTile,
          typename Index, typename N, typename F,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
void transformParticles (DstTile1& dst1, DstTile2& dst2, const SrcTile& src,
                         Index src_start, Index dst1_start, Index dst2_start, N n, F const& f) noexcept
{
    const auto src_data  = src.getConstParticleTileData();
          auto dst1_data = dst1.getParticleTileData();
          auto dst2_data = dst2.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        f(dst1_data, dst2_data, src_data, src_start+i, dst1_start+i, dst2_start+i);
    });

    Gpu::streamSynchronize();
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename N,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
Index filterParticles (DstTile& dst, const SrcTile& src, const Index* mask) noexcept
{
    return filterParticles(dst, src, mask, 0, 0, src.numParticles());
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 *  This version conditionally copies n particles starting at index src_start, writing
 *  the result starting at dst_start.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param n the number of particles to apply the operation to
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename N,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
Index filterParticles (DstTile& dst, const SrcTile& src, const Index* mask,
                       Index src_start, Index dst_start, N n) noexcept
{
    Gpu::DeviceVector<Index> offsets(n);
    Gpu::exclusive_scan(mask, mask+n, offsets.begin());

    Index last_mask=0, last_offset=0;
    Gpu::copyAsync(Gpu::deviceToHost, mask+n-1, mask + n, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+n-1, offsets.data()+n, &last_offset);

    auto* p_offsets = offsets.dataPtr();

    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( n, i,
    {
        if (mask[i]) { copyParticle(dst_data, src_data, src_start+i, dst_start+p_offsets[i]); }
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 *
 */
template <typename DstTile, typename SrcTile, typename Pred,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<Pred>>,int> foo = 0>
int filterParticles (DstTile& dst, const SrcTile& src, Pred&& p) noexcept
{
    return filterParticles(dst, src, std::forward<Pred>(p), 0, 0, src.numParticles());
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 *  This version conditionally copies n particles starting at index src_start, writing
 *  the result starting at dst_start.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 * \param n the number of particles to apply the operation to
 *
 */
template <typename DstTile, typename SrcTile, typename Pred, typename Index, typename N,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<Pred>>,Index> nvccfoo = 0>
Index filterParticles (DstTile& dst, const SrcTile& src, Pred const& p,
                       Index src_start, Index dst_start, N n) noexcept
{
    Gpu::DeviceVector<Index> mask(n);

    auto* p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();

    amrex::ParallelForRNG(n,
    [p, p_mask, src_data, src_start] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
    {
        amrex::ignore_unused(p, p_mask, src_data, src_start, engine);
        if constexpr (IsCallable<Pred,decltype(src_data),Index,RandomEngine>::value) {
            p_mask[i] = p(src_data, src_start+i, engine);
        } else {
            p_mask[i] = p(src_data, src_start+i);
        }
    });
    return filterParticles(dst, src, mask.dataPtr(), src_start, dst_start, n);
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename F,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
Index filterAndTransformParticles (DstTile& dst, const SrcTile& src, Index* mask, F const& f,
                                   Index src_start, Index dst_start) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask=0, last_offset=0;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    auto const* p_offsets = offsets.dataPtr();

    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i]) {
            f(dst_data, src_data, src_start+i,
              dst_start+p_offsets[src_start+i]);
        }
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile, typename SrcTile, typename Index, typename F,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
Index filterAndTransformParticles (DstTile& dst, const SrcTile& src, Index* mask, F&& f) noexcept
{
    return filterAndTransformParticles(dst, src, mask, std::forward<F>(f), 0, 0);
}

/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 * \tparam F the transform function type
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile, typename SrcTile, typename Pred, typename F,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<Pred>>,int> foo = 0>
int filterAndTransformParticles (DstTile& dst, const SrcTile& src, Pred&& p, F&& f) noexcept
{
    return filterAndTransformParticles(dst, src, std::forward<Pred>(p), std::forward<F>(f), 0, 0);
}

/**
 * \brief Conditionally copy particles from src to dst1 and dst2 based on the value of mask.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile1 the dst1 particle tile type
 * \tparam DstTile2 the dst2 particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Index the index type, e.g. unsigned int
 * \tparam F the transform function type
 *
 * \param dst1 the first destination tile
 * \param dst2 the second destination tile
 * \param src the source tile
 * \param mask pointer to the mask - 1 means copy, 0 means don't copy
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile1, typename DstTile2, typename SrcTile, typename Index, typename F,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
Index filterAndTransformParticles (DstTile1& dst1, DstTile2& dst2,
                                   const SrcTile& src, Index* mask, F const& f) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<Index> offsets(np);
    Gpu::exclusive_scan(mask, mask+np, offsets.begin());

    Index last_mask=0, last_offset=0;
    Gpu::copyAsync(Gpu::deviceToHost, mask+np-1, mask + np, &last_mask);
    Gpu::copyAsync(Gpu::deviceToHost, offsets.data()+np-1, offsets.data()+np, &last_offset);

    auto* p_offsets = offsets.dataPtr();

    const auto src_data  = src.getConstParticleTileData();
          auto dst_data1 = dst1.getParticleTileData();
          auto dst_data2 = dst2.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        if (mask[i]) { f(dst_data1, dst_data2, src_data, i, p_offsets[i], p_offsets[i]); }
    });

    Gpu::streamSynchronize();
    return last_mask + last_offset;
}

/**
 * \brief Conditionally copy particles from src to dst1 and dst2 based on a predicate.
 * A transformation will also be applied to the particles on copy.
 *
 * \tparam DstTile1 the dst1 particle tile type
 * \tparam DstTile2 the dst2 particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 * \tparam F the transform function type
 *
 * \param dst1 the first destination tile
 * \param dst2 the second destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param f defines the transformation that will be applied to the particles on copy
 *
 */
template <typename DstTile1, typename DstTile2, typename SrcTile, typename Pred, typename F,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<Pred>>, int> foo = 0>
int filterAndTransformParticles (DstTile1& dst1, DstTile2& dst2, const SrcTile& src,
                                 Pred const& p, F&& f) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<int> mask(np);

    auto* p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();

    amrex::ParallelForRNG(np,
    [p, p_mask, src_data] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
    {
        amrex::ignore_unused(p, p_mask, src_data, engine);
        if constexpr (IsCallable<Pred,decltype(src_data),int,RandomEngine>::value) {
            p_mask[i] = p(src_data, i, engine);
        } else {
            p_mask[i] = p(src_data, i);
        }
    });
    return filterAndTransformParticles(dst1, dst2, src, mask.dataPtr(), std::forward<F>(f));
}


/**
 * \brief Conditionally copy particles from src to dst based on a predicate.
 *  This version conditionally copies n particles starting at index src_start, writing
 *  the result starting at dst_start.
 *
 * \tparam DstTile the dst particle tile type
 * \tparam SrcTile the src particle tile type
 * \tparam Pred a function object
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param p predicate function - particles will be copied if p returns true
 * \param src_start the offset at which to start reading particles from src
 * \param dst_start the offset at which to start writing particles to dst
 *
 */
template <typename DstTile, typename SrcTile, typename Pred, typename F, typename Index,
          std::enable_if_t<!std::is_pointer_v<std::decay_t<Pred>>,Index> nvccfoo = 0>
Index filterAndTransformParticles (DstTile& dst, const SrcTile& src, Pred const& p, F&& f,
                                   Index src_start, Index dst_start) noexcept
{
    auto np = src.numParticles();
    Gpu::DeviceVector<Index> mask(np);

    auto* p_mask = mask.dataPtr();
    const auto src_data = src.getConstParticleTileData();

    amrex::ParallelForRNG(np,
    [p, p_mask, src_data, src_start] AMREX_GPU_DEVICE (int i, amrex::RandomEngine const& engine) noexcept
    {
        amrex::ignore_unused(p, p_mask, src_data, src_start, engine);
        if constexpr (IsCallable<Pred,decltype(src_data),Index,RandomEngine>::value) {
            p_mask[i] = p(src_data, src_start+i, engine);
        } else {
            p_mask[i] = p(src_data, src_start+i);
        }
    });
    return filterAndTransformParticles(dst, src, mask.dataPtr(), std::forward<F>(f), src_start, dst_start);
}


/**
 * \brief Gather particles copies particles into contiguous order from an
 * arbitrary order. Specifically, the particle at the index inds[i] in src
 * will be copied to the index i in dst.
 *
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. Long
 * \tparam Index the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename Index,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
void gatherParticles (PTile& dst, const PTile& src, N np, const Index* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        copyParticle(dst_data, src_data, inds[i], i);
    });

    Gpu::streamSynchronize();
}

/**
 * \brief Scatter particles copies particles from contiguous order into an
 * arbitrary order. Specifically, the particle at the index i in src
 * will be copied to the index inds[i] in dst.
 *
 * \tparam PTile the particle tile type
 * \tparam N the size type, e.g. Long
 * \tparam Index the index type, e.g. unsigned int
 *
 * \param dst the destination tile
 * \param src the source tile
 * \param np the number of particles
 * \param inds pointer to the permutation array
 *
 */
template <typename PTile, typename N, typename Index,
          std::enable_if_t<std::is_integral_v<Index>, int> foo = 0>
void scatterParticles (PTile& dst, const PTile& src, N np, const Index* inds)
{
    const auto src_data = src.getConstParticleTileData();
          auto dst_data = dst.getParticleTileData();

    AMREX_HOST_DEVICE_FOR_1D( np, i,
    {
        copyParticle(dst_data, src_data, i, inds[i]);
    });

    Gpu::streamSynchronize();
}

}

#endif // include guard
