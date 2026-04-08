#ifndef AMREX_GPU_TENSOR_COPY_H_
#define AMREX_GPU_TENSOR_COPY_H_
#include <AMReX_Config.H>

#include <AMReX_Box.H>
#include <AMReX_GpuDevice.H>

#include <cstddef>
#include <cstdint>

#if defined(AMREX_USE_CUDA) && defined(__CUDACC__) && defined(__has_include)
#  if __has_include(<cub/block/block_load_to_shared.cuh>) && __has_include(<cuda/std/span>)
#    define AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED 1
#  endif
#endif

#ifndef AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
#  define AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED 0
#endif

namespace amrex::Gpu::detail {

struct TensorCopyCaps {
    bool available = false;
    bool requires_alignment_check = true;
    bool supports_multidim = false;
    bool supports_shared_padding = false;
};

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
constexpr std::size_t tensorCopyAlignment () noexcept
{
    return 16u;
}

[[nodiscard]] AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool isTensorCopyAligned (void const* p, std::size_t alignment = tensorCopyAlignment()) noexcept
{
    return (reinterpret_cast<std::uintptr_t>(p) % alignment) == 0u;
}

[[nodiscard]] inline TensorCopyCaps queryTensorCopyCaps () noexcept
{
    TensorCopyCaps caps;
#if AMREX_USE_CUDA_BLOCK_LOAD_TO_SHARED
    caps.available = (Gpu::Device::devicePropMajor() >= 9);
    caps.supports_multidim = true;
    caps.supports_shared_padding = true;
#endif
    return caps;
}

[[nodiscard]] inline bool canUseTensorCopy (void const* p,
                                            Box const& src_box,
                                            IntVect const& tile,
                                            IntVect const& halo,
                                            std::size_t element_size) noexcept
{
    amrex::ignore_unused(src_box);

    TensorCopyCaps const caps = queryTensorCopyCaps();
    if (!caps.available || (element_size == 0u)) {
        return false;
    }

    IntVect const staged = tile + halo + halo;
    if (staged.min() <= 0) {
        return false;
    }

    Dim3 const staged_d = staged.dim3(1);
    auto const nelems = static_cast<std::size_t>(staged_d.x)
        * static_cast<std::size_t>(staged_d.y)
        * static_cast<std::size_t>(staged_d.z);

#ifdef AMREX_USE_GPU
    return isTensorCopyAligned(p)
        && nelems * element_size <= Gpu::Device::sharedMemPerBlock();
#else
    amrex::ignore_unused(nelems);
    return false;
#endif
}

} // namespace amrex::Gpu::detail

namespace amrex::Gpu {

[[nodiscard]] inline bool supportsTensorCopy () noexcept
{
    return detail::queryTensorCopyCaps().available;
}

} // namespace amrex::Gpu

#endif
