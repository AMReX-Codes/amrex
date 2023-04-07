#ifndef AMREX_GPU_CONTAINERS_H_
#define AMREX_GPU_CONTAINERS_H_
#include <AMReX_Config.H>

#include <AMReX_Vector.H>
#include <AMReX_PODVector.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_Scan.H>
#include <type_traits>

#include <numeric>
#include <iterator>

namespace amrex::Gpu {

#ifdef AMREX_USE_GPU

    /**
     * \brief A PODVector that uses the standard memory Arena.
     * Note that the memory might or might not be managed depending
     * on the amrex.the_arena_is_managed ParmParse parameter.
     */
    template <class T>
    using DeviceVector = PODVector<T, ArenaAllocator<T> >;

    /**
     * \brief A PODVector that uses the non-managed device memory arena.
     *
     */
    template <class T>
    using NonManagedDeviceVector = PODVector<T, DeviceArenaAllocator<T> >;

    /**
     * \brief A PODVector that uses the managed memory arena.
     *
     */
    template <class T>
    using ManagedVector = PODVector<T, ManagedArenaAllocator<T> >;

    /**
     * \brief A PODVector that uses the pinned memory arena.
     *
     */
    template <class T>
    using PinnedVector = PODVector<T, PinnedArenaAllocator<T> >;

    /**
     * \brief A PODVector that uses the async memory arena.
     *        Maybe useful for temporary vectors inside MFIters
     *        that are accessed on the device.
     *
     */
    template <class T>
    using AsyncVector = PODVector<T, AsyncArenaAllocator<T> >;

    /**
     * \brief A PODVector that uses pinned host memory. Same as PinnedVector.
     * For a vector class that uses std::allocator by default, see amrex::Vector.
     *
     */
    template <class T>
    using HostVector = PinnedVector<T>;

    /**
     * \brief The behavior of PolymorphicVector changes depending on
     * the amrex.use_gpu_aware_mpi runtime flag. If the flag is true,
     * this vector will use device memory. If it is false, this Vector
     * will use pinned memory.
     *
     */
    template <class T>
    using PolymorphicVector = PODVector<T, PolymorphicAllocator<T> >;

    /**
     * \brief This is identical to ManagedVector<T>. The ManagedDeviceVector
     * form is deprecated and will be removed in a future release.
     *
     */
    template <class T>
    using ManagedDeviceVector = PODVector<T, ManagedArenaAllocator<T> >;

#else
    //! When Cuda is off, all these containers revert to amrex::Gpu::PODVector.
    template <class T>
    using DeviceVector = PODVector<T>;

    template <class T>
    using HostVector = PODVector<T>;

    template <class T>
    using NonManagedVector = PODVector<T>;

    template <class T>
    using ManagedVector = PODVector<T>;

    template <class T>
    using ManagedDeviceVector = PODVector<T>;

    template <class T>
    using PinnedVector = PODVector<T>;

    template <class T>
    using AsyncVector = PODVector<T>;

    template <class T>
    using PolymorphicVector = PODVector<T>;
#endif

    struct HostToDevice {};
    struct DeviceToHost {};
    struct DeviceToDevice {};
    static constexpr HostToDevice   hostToDevice{};
    static constexpr DeviceToHost   deviceToHost{};
    static constexpr DeviceToDevice deviceToDevice{};

    /**
     * \brief A host-to-device copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is blocking - CPU execution will halt until the copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param HostToDevice A tag indicating that the copy is from the host to the device
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copy(Gpu::hostToDevice, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copy (HostToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        htod_memcpy(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief A device-to-host copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is blocking - CPU execution will halt until the copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param DeviceToHost A tag indicating that the copy is from the device to the host
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copy(Gpu::deviceToHost, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copy (DeviceToHost, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        dtoh_memcpy(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief A device-to-device copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is blocking - CPU execution will halt until the copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param DeviceToDevice A tag indicating that the copy is from the device to the device
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copy(Gpu::deviceToDevice, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copy (DeviceToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        dtod_memcpy(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief A host-to-device copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is asynchronous - CPU execution will continue, whether or not the
     * copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param HostToDevice A tag indicating that the copy is from the host to the device
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copyAsync(Gpu::hostToDevice, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copyAsync (HostToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        htod_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief A device-to-host copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is asynchronous - CPU execution will continue, whether or not the
     * copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param DeviceToHost A tag indicating that the copy is from the device to the host
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copyAsync(Gpu::deviceToHost, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copyAsync (DeviceToHost, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        dtoh_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief A device-to-device copy routine. Note this is just a wrapper around
     * memcpy, so it assumes contiguous storage. The amrex-provided containers
     * like Gpu::HostVector, Gpu::DeviceVector, etc. meet this requirement.
     *
     * This version is asynchronous - CPU execution will continue, whether or not the
     * copy is finished.
     *
     * \tparam InIter  The input iterator type
     * \tparam OutIter The output iterator type
     * \param DeviceToDevice A tag indicating that the copy is from the device to the device
     * \param begin Where in the input to start reading
     * \param end Where in the input to stop reading
     * \param result Where in the output to start writing
     *
     * Example usage:
     *
     *    Gpu::copyAsync(Gpu::deviceToDevice, a.begin(), a.end(), b.begin());
     */
    template<class InIter, class OutIter>
    void copyAsync (DeviceToDevice, InIter begin, InIter end, OutIter result) noexcept
    {
        using value_type = typename std::iterator_traits<InIter>::value_type;

        using out_value_type = typename std::iterator_traits<OutIter>::value_type;
        static_assert(std::is_same_v<value_type, out_value_type>);
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;
#ifdef AMREX_USE_GPU
        dtod_memcpy_async(&(*result), &(*begin), size*sizeof(value_type));
#else
        std::memcpy(&(*result), &(*begin), size*sizeof(value_type));
#endif
    }

    /**
     * \brief Migrate elements of a container from device to host.
     * This is a no-op for host-only code.
     *
     * This version is blocking - CPU execution will halt until the migration is finished.
     */
    template<class Iter>
    void prefetchToHost (Iter begin, Iter end) noexcept
    {
        using value_type = typename std::iterator_traits<Iter>::value_type;
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;

#ifdef AMREX_USE_GPU
        // Currently only implemented for CUDA.
#if defined(AMREX_USE_CUDA) && !defined(_WIN32)
        if (Gpu::Device::devicePropMajor() >= 6) {
            AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(&(*begin),
                                                      size*sizeof(value_type),
                                                      cudaCpuDeviceId,
                                                      Gpu::gpuStream()));
        }
#endif
#endif

        Gpu::streamSynchronize();
    }

    /**
     * \brief Migrate elements of a container from host to device.
     * This is a no-op for host-only code.
     *
     * This version is blocking - CPU execution will halt until the migration is finished.
     */
    template<class Iter>
    void prefetchToDevice (Iter begin, Iter end) noexcept
    {
        using value_type = typename std::iterator_traits<Iter>::value_type;
        static_assert(std::is_trivially_copyable<value_type>(),
                      "Can only copy trivially copyable types");

        auto size = std::distance(begin, end);
        if (size == 0) return;

#ifdef AMREX_USE_GPU
        // Currently only implemented for CUDA.
#if defined(AMREX_USE_CUDA) && !defined(_WIN32)
        if (Gpu::Device::devicePropMajor() >= 6) {
            AMREX_CUDA_SAFE_CALL(cudaMemPrefetchAsync(&(*begin),
                                                      size*sizeof(value_type),
                                                      Gpu::Device::deviceId(),
                                                      Gpu::gpuStream()));
        }
#endif
#endif

        Gpu::streamSynchronize();
    }

    /**
     * \brief Fill the elements in the given range using the given
     * calllable.
     *
     *  This function is asynchronous for GPU builds.
     *
     * \tparam IT the iterator type
     * \tparam F  the callable type
     *
     * \param first the inclusive first in the range [first, last)
     * \param last  the exclusive last in the range [first, last)
     * \param f     the callable with the function signature of void(T&, Long),
     *              where T is the element type and the Long parameter is the
     *              index for the element to be filled.
     */
    template <typename IT, typename F,
              typename T = typename std::iterator_traits<IT>::value_type,
              std::enable_if_t<(sizeof(T) <= 36*8) &&  // so there is enough shared memory
                               std::is_trivially_copyable_v<T> &&
                               amrex::IsCallable<F, T&, Long>::value,
                               int> FOO = 0>
    void fillAsync (IT first, IT last, F&& f) noexcept
    {
        auto N = static_cast<Long>(std::distance(first, last));
        if (N <= 0) return;
        auto p = &(*first);
#ifndef AMREX_USE_GPU
        for (Long i = 0; i < N; ++i) {
            f(p[i], i);
        }
#else
        // No need to use shared memory if the type is small.
        // May not have enough shared memory if the type is too big.
        // Cannot use shared memory, if the type is not trivially copable.
        if constexpr ((sizeof(T) <= 8)
                      || (sizeof(T) > 36*8)
                      || ! std::is_trivially_copyable<T>()) {
            amrex::ParallelFor(N, [=] AMREX_GPU_DEVICE (Long i) noexcept
            {
                f(p[i], i);
            });
        } else {
            static_assert(sizeof(T) % sizeof(unsigned int) == 0);
            using U = std::conditional_t<sizeof(T) % sizeof(unsigned long long) == 0,
                                         unsigned long long, unsigned int>;
            constexpr Long nU = sizeof(T) / sizeof(U);
            auto pu = reinterpret_cast<U*>(p);
            int nthreads_per_block = (sizeof(T) <= 64) ? 256 : 128;
            int nblocks = static_cast<int>((N+nthreads_per_block-1)/nthreads_per_block);
            std::size_t shared_mem_bytes = nthreads_per_block * sizeof(T);
#ifdef AMREX_USE_SYCL
            amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
            [=] AMREX_GPU_DEVICE (Gpu::Handler const& handler) noexcept
            {
                Long i = handler.globalIdx();
                Long blockDimx = handler.blockDim();
                Long threadIdxx = handler.threadIdx();
                Long blockIdxx = handler.blockIdx();
                auto const shared_U = (U*)handler.sharedMemory();
                auto const shared_T = (T*)shared_U;
                if (i < N) {
                    auto ga = new(shared_T+threadIdxx) T;
                    f(*ga, i);
                }
                handler.sharedBarrier();
                for (Long m = threadIdxx,
                         mend = nU * amrex::min(blockDimx, N-blockDimx*blockIdxx);
                     m < mend; m += blockDimx) {
                    pu[blockDimx*blockIdxx*nU+m] = shared_U[m];
                }
            });
#else
            amrex::launch(nblocks, nthreads_per_block, shared_mem_bytes, Gpu::gpuStream(),
                          [=] AMREX_GPU_DEVICE () noexcept
            {
                Long blockDimx = blockDim.x;
                Long threadIdxx = threadIdx.x;
                Long blockIdxx = blockIdx.x;
                Long i = blockDimx*blockIdxx + threadIdxx;
                Gpu::SharedMemory<U> gsm;
                auto const shared_U = gsm.dataPtr();
                auto const shared_T = (T*)shared_U;
                if (i < N) {
                    auto ga = new(shared_T+threadIdxx) T;
                    f(*ga, i);
                }
                __syncthreads();
                for (Long m = threadIdxx,
                         mend = nU * amrex::min(blockDimx, N-blockDimx*blockIdxx);
                     m < mend; m += blockDimx) {
                    pu[blockDimx*blockIdxx*nU+m] = shared_U[m];
                }
            });
#endif
        }
#endif
    }

}

#endif
