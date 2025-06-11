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
     * Given the diverse interfaces of FFT libraries we use, this constructo
     * has a number of optional arguments.
     *
     * The user can provide the data pointers to the constructor. They are
     * only needed by FFTW because its plan creation requires the input and
     * output arrays. If they are null, we will delay the plan creation for
     * FFTW until the forward or backward function is called.
     *
     * The cache_plan option is only used when we use cufft, rocfft and
     * onemkl, but not FFTW.
     *
     * \param fft_size The forward domain size (i.e., the domain of the real data)
     * \param p_fwd Forward domain data pointer (optional)
     * \param p_bwd Backward domain data pointer (optional)
     * \param cache_plan Try to cache the plan or not (optionl)
     */
    explicit LocalR2C (IntVectND<M> const& fft_size,
                       T* p_fwd = nullptr,
                       GpuComplex<T>* p_bwd = nullptr,
#ifdef AMREX_USE_GPU
                       bool cache_plan = true);
#else
                       bool cache_plan = false);
#endif

    ~LocalR2C ();

    LocalR2C () = default;
    LocalR2C (LocalR2C &&) noexcept;
    LocalR2C& operator= (LocalR2C &&) noexcept;

    LocalR2C (LocalR2C const&) = delete;
    LocalR2C& operator= (LocalR2C const&) = delete;

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

    void clear ();

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

    T* m_p_fwd = nullptr;
    GpuComplex<T>* m_p_bwd = nullptr;

#if defined(AMREX_USE_SYCL)
    gpuStream_t m_gpu_stream{};
#endif

    IntVectND<M> m_real_size;
    IntVectND<M> m_spectral_size;

    bool m_cache_plan = false;
};

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>::LocalR2C (IntVectND<M> const& fft_size, T* p_fwd,
                           GpuComplex<T>* p_bwd, bool cache_plan)
    : m_p_fwd(p_fwd),
      m_p_bwd(p_bwd),
      m_real_size(fft_size),
      m_spectral_size(fft_size)
#if defined(AMREX_USE_GPU)
      , m_cache_plan(cache_plan)
#endif
{
#if !defined(AMREX_USE_GPU)
    amrex::ignore_unused(cache_plan);
#endif

    BL_PROFILE("FFT::LocalR2C");
    m_spectral_size[0] = m_real_size[0]/2 + 1;

#if defined(AMREX_USE_SYCL)

    auto current_stream = Gpu::gpuStream();
    Gpu::Device::resetStreamIndex();
    m_gpu_stream = Gpu::gpuStream();

#endif

    auto* pf = (void*)m_p_fwd;
    auto* pb = (void*)m_p_bwd;

#ifdef AMREX_USE_SYCL
    m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, pf, pb, m_cache_plan);
    m_fft_bwd = m_fft_fwd;
#else
    if constexpr (D == Direction::both || D == Direction::forward) {
        m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, pf, pb, m_cache_plan);
    }
    if constexpr (D == Direction::both || D == Direction::backward) {
        m_fft_bwd.template init_r2c<Direction::backward,M>(m_real_size, pf, pb, m_cache_plan);
    }
#endif

#if defined(AMREX_USE_SYCL)
    Gpu::Device::setStream(current_stream);
#endif
}

template <typename T, FFT::Direction D, int M>
void LocalR2C<T,D,M>::clear ()
{
    if (!m_cache_plan) {
        if (m_fft_bwd.plan != m_fft_fwd.plan) {
            m_fft_bwd.destroy();
        }
        m_fft_fwd.destroy();
    }

    m_fft_fwd = Plan<T>{};
    m_fft_bwd = Plan<T>{};
}

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>::~LocalR2C ()
{
    static_assert(M >= 1 && M <= 3);
    clear();
}

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>::LocalR2C (LocalR2C && rhs) noexcept
    : m_p_fwd(rhs.m_p_fwd),
      m_p_bwd(rhs.m_p_bwd),
      m_fft_fwd(rhs.m_fft_fwd),
      m_fft_bwd(rhs.m_fft_bwd),
#if defined(AMREX_USE_SYCL)
      m_gpu_stream(rhs.m_gpu_stream),
#endif
      m_real_size(rhs.m_real_size),
      m_spectral_size(rhs.m_spectral_size),
      m_cache_plan(rhs.m_cache_plan)
{
    rhs.m_cache_plan = true; // So that plans in rhs are not destroyed.
}

template <typename T, FFT::Direction D, int M>
LocalR2C<T,D,M>& LocalR2C<T,D,M>::operator= (LocalR2C && rhs) noexcept
{
    if (this == &rhs) { return *this; }

    this->clear();

    m_p_fwd = rhs.m_p_fwd;
    m_p_bwd = rhs.m_p_bwd;
    m_fft_fwd = rhs.m_fft_fwd;
    m_fft_bwd = rhs.m_fft_bwd;
#if defined(AMREX_USE_SYCL)
    m_gpu_stream = rhs.m_gpu_stream;
#endif
    m_real_size = rhs.m_real_size;
    m_spectral_size = rhs.m_spectral_size;
    m_cache_plan = rhs.m_cache_plan;

    rhs.m_cache_plan = true; // So that plans in rhs are not destroyed.

    return *this;
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

    if (((T*)indata != m_p_fwd) || (outdata != m_p_bwd)) {
        m_p_fwd = (T*)indata;
        m_p_bwd = outdata;
        auto* pf = (void*)m_p_fwd;
        auto* pb = (void*)m_p_bwd;
        m_fft_fwd.destroy();
        m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, pf, pb, false);
        if constexpr (D == Direction::both) {
            m_fft_bwd.destroy();
            m_fft_bwd.template init_r2c<Direction::backward,M>(m_real_size, pf, pb, false);
        }
    }

#endif

    m_fft_fwd.template compute_r2c<Direction::forward>();

#if defined(AMREX_USE_SYCL)
    if (current_stream != m_gpu_stream) {
        Gpu::Device::setStream(current_stream);
    }
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

    if (((GpuComplex<T>*)indata != m_p_bwd) || (outdata != m_p_fwd)) {
        m_p_fwd = outdata;
        m_p_bwd = (GpuComplex<T>*)indata;
        auto* pf = (void*)m_p_fwd;
        auto* pb = (void*)m_p_bwd;
        m_fft_bwd.destroy();
        m_fft_bwd.template init_r2c<Direction::backward,M>(m_real_size, pf, pb, false);
        if constexpr (D == Direction::both) {
            m_fft_fwd.destroy();
            m_fft_fwd.template init_r2c<Direction::forward,M>(m_real_size, pf, pb, false);
        }
    }

#endif

    m_fft_bwd.template compute_r2c<Direction::backward>();

#if defined(AMREX_USE_SYCL)
    if (current_stream != m_gpu_stream) {
        Gpu::Device::setStream(current_stream);
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
