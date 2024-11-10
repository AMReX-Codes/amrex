#ifndef AMREX_FFT_LOCAL_R2C_H_
#define AMREX_FFT_LOCAL_R2C_H_
#include <AMReX_Config.H>

#include <AMReX_Arena.H>
#include <AMReX_FFT_Helper.H>

namespace amrex::FFT
{

/**
 * \brief Local Discrete Fourier Transform
 *
 * This class supports Fourier transforms between real and complex data. The
 * name R2C indicates that the forward transform converts real data to
 * complex data, while the backward transform converts complex data to real
 * data. It should be noted that both directions of transformation are
 * supported, not just from real to complex. The scaling follows the FFTW
 * convention, where applying the forward transform followed by the backward
 * transform scales the original data by the size of the input array.
 *
 * For more details, we refer the users to
 * https://amrex-codes.github.io/amrex/docs_html/FFT_Chapter.html.
 */
template <typename T, FFT::Direction D = FFT::Direction::both,
          int M = AMREX_SPACEDIM>
class LocalR2C
{
public:
    /**
     * \brief Constructor
     *
     * \param domain the forward domain (i.e., the domain of the real data)
     * \param info optional information
     */
    explicit LocalR2C (IntVectND<M> const& fft_size, bool cache_plan = true);

    ~LocalR2C ();

    LocalR2C (LocalR2C const&) = delete;
    LocalR2C (LocalR2C &&) = delete;
    LocalR2C& operator= (LocalR2C const&) = delete;
    LocalR2C& operator= (LocalR2C &&) = delete;

    /**
     * \brief Forward transform
     *
     * This function is not available when this class template is
     * instantiated for backward-only transform. For GPUs, this function is
     * synchronous on the host.
     *
     * \param indata input data
     * \param outdata output data
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::forward ||
                                                DIR == Direction::both, int> = 0>
    void forward (T const* indata, GpuComplex<T>* outdata);

    /**
     * \brief Backward transform
     *
     * This function is not available when this class template is
     * instantiated for forward-only transform. For GPUs, this function is
     * synchronous on the host.
     *
     * \param indata input data
     * \param outdata output data
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::backward ||
                                                DIR == Direction::both, int> = 0>
    void backward (GpuComplex<T> const* indata, T* outdata);

    //! Scaling factor. If the data goes through forward and then backward,
    //! the result multiplied by the scaling factor is equal to the original
    //! data.
    [[nodiscard]] T scalingFactor () const;

    //! Spectral domain size
    [[nodiscard]] IntVectND<M> const& spectralSize () const {
        return m_spectral_size;
    }

private:

    Plan<T> m_fft_fwd;
    Plan<T> m_fft_bwd;

    void* m_buffer = nullptr; // for fftw

#if defined(AMREX_USE_SYCL)
    gpuStream_t m_gpu_stream;
#endif

    IntVectND<M> m_real_size;
    IntVectND<M> m_spectral_size;

    bool m_cache_plan;
};

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>::LocalR2C (IntVectND<M> const& fft_size, bool cache_plan)
    : m_real_size(fft_size),
      m_spectral_size(fft_size),
      m_cache_plan(cache_plan)
{
    BL_PROFILE("FFT::LocalR2C");
    m_spectral_size[0] = m_real_size[0]/2 + 1;

#if defined(AMREX_USE_SYCL)

    auto current_stream = Gpu::gpuStream();
    Gpu::Device::resetStreamIndex();
    m_gpu_stream = Gpu::gpuStream();

#elif !defined(AMREX_USE_GPU)

    Long total_size = 1;
    for (auto s : m_spectral_size) {
        total_size *= s;
    }
    m_buffer = The_Arena()->alloc(sizeof(GpuComplex<T>)*total_size);

#endif

#ifdef AMREX_USE_SYCL
    m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, m_buffer, cache_plan);
    m_fft_bwd = m_fft_fwd;
#else
    if constexpr (D == Direction::both || D == Direction::forward) {
        m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, m_buffer, cache_plan);
    }
    if constexpr (D == Direction::both || D == Direction::backward) {
        m_fft_bwd.template init_r2c<Direction::backward,M>(m_real_size, m_buffer, cache_plan);
    }
#endif

#if defined(AMREX_USE_SYCL)
    Gpu::Device::setStream(current_stream);
#endif
}

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>::~LocalR2C ()
{
    static_assert(M >= 1 && M <= 3);

    if (m_buffer) { The_Arena()->free(m_buffer); }

    if (!m_cache_plan) {
        if (m_fft_bwd.plan != m_fft_fwd.plan) {
            m_fft_bwd.destroy();
        }
        m_fft_fwd.destroy();
    }
}

template <typename T, FFT::Direction D, int M>
template <Direction DIR, std::enable_if_t<DIR == Direction::forward ||
                                          DIR == Direction::both, int> >
void LocalR2C<T,D,M>::forward (T const* indata, GpuComplex<T>* outdata)
{
    BL_PROFILE("FFT::LocalR2C::forward");

#if defined(AMREX_USE_GPU)

    m_fft_fwd.set_ptrs((void*)indata, (void*)outdata);

#if defined(AMREX_USE_SYCL)
    auto current_stream = Gpu::gpuStream();
    if (current_stream != m_gpu_stream) {
        Gpu::streamSynchronize();
        Gpu::Device::setStream(m_gpu_stream);
    }
#endif

#else /* FFTW */

    int ny = 1;
    for (int idim = 1; idim < M; ++idim) {
        ny *= m_real_size[idim];
    }
    auto* pdst = (char*)m_buffer;
    auto* psrc = (char*)indata;
    for (int j = 0; j < ny; ++j) {
        std::memcpy(pdst, psrc, sizeof(T) * m_real_size[0]);
        pdst += sizeof(GpuComplex<T>) * m_spectral_size[0];
        psrc += sizeof(           T ) *     m_real_size[0];
    }

#endif

    m_fft_fwd.template compute_r2c<Direction::forward>();

#if defined(AMREX_USE_SYCL)

    if (current_stream != m_gpu_stream) {
        Gpu::Device::setStream(current_stream);
    }

#elif !defined(AMREX_USE_GPU)

    std::size_t nbytes = sizeof(GpuComplex<T>);
    for (auto s : m_spectral_size) { nbytes *= s; }
    std::memcpy(outdata, m_buffer, nbytes);

#endif
}

template <typename T, FFT::Direction D, int M>
template <Direction DIR, std::enable_if_t<DIR == Direction::backward ||
                                          DIR == Direction::both, int> >
void LocalR2C<T,D,M>::backward (GpuComplex<T> const* indata, T* outdata)
{
    BL_PROFILE("FFT::LocalR2C::backward");

#if defined(AMREX_USE_GPU)

    m_fft_bwd.set_ptrs((void*)outdata, (void*)indata);

#if defined(AMREX_USE_SYCL)
    auto current_stream = Gpu::gpuStream();
    if (current_stream != m_gpu_stream) {
        Gpu::streamSynchronize();
        Gpu::Device::setStream(m_gpu_stream);
    }
#endif

#else /* FFTW */

    std::size_t nbytes = sizeof(GpuComplex<T>);
    for (auto s : m_spectral_size) { nbytes *= s; }
    std::memcpy(m_buffer, indata, nbytes);

#endif

    m_fft_bwd.template compute_r2c<Direction::backward>();

#if defined(AMREX_USE_SYCL)

    if (current_stream != m_gpu_stream) {
        Gpu::Device::setStream(current_stream);
    }

#elif !defined(AMREX_USE_GPU)

    int ny = 1;
    for (int idim = 1; idim < M; ++idim) {
        ny *= m_real_size[idim];
    }
    auto* pdst = (char*)outdata;
    auto* psrc = (char*)m_buffer;
    for (int j = 0; j < ny; ++j) {
        std::memcpy(pdst, psrc, sizeof(T) * m_real_size[0]);
        pdst += sizeof(           T ) *     m_real_size[0];
        psrc += sizeof(GpuComplex<T>) * m_spectral_size[0];
    }

#endif
}

template <typename T, FFT::Direction D, int M>
T LocalR2C<T,D,M>::scalingFactor () const
{
    T r = 1;
    for (auto s : m_real_size) {
        r *= T(s);
    }
    return T(1)/r;
}

}

#endif
