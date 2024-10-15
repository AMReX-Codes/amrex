#ifndef AMREX_FFT_H_
#define AMREX_FFT_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <numeric>
#include <tuple>
#include <utility>

#if defined(AMREX_USE_CUDA)
#  include <cufft.h>
#  include <cuComplex.h>
#elif defined(AMREX_USE_HIP)
#  if __has_include(<rocfft/rocfft.h>)  // ROCm 5.3+
#    include <rocfft/rocfft.h>
#  else
#    include <rocfft.h>
#  endif
#  include <hip/hip_complex.h>
#elif defined(AMREX_USE_SYCL)
#  include <oneapi/mkl/dfti.hpp>
#else
#  include <fftw3.h>
#endif

namespace amrex::FFT
{

enum struct Scaling { full, symmetric, none };
enum struct Direction { forward, backward };

template <typename T = Real>
class R2C
{
public:
    using MF = std::conditional_t<std::is_same_v<T,Real>,
                                  MultiFab, FabArray<BaseFab<T> > >;
    using cMF = FabArray<BaseFab<GpuComplex<T> > >;

    R2C (Box const& domain);

    ~R2C ();

    R2C (R2C const&) = delete;
    R2C (R2C &&) = delete;
    R2C& operator= (R2C const&) = delete;
    R2C& operator= (R2C &&) = delete;

    template <typename F>
    void forwardThenBackward (MF const& inmf, MF& outmf, F const& post_forward)
    {
        this->forward_doit(inmf);
        this->post_forward_doit(post_forward);
        this->backward_doit(outmf);
    }

    struct Swap01
    {
        [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 operator() (Dim3 i) const noexcept
        {
            return {i.y, i.x, i.z};
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 Inverse (Dim3 i) const noexcept
        {
            return {i.y, i.x, i.z};
        }

        [[nodiscard]] IndexType operator() (IndexType it) const noexcept
        {
            return it;
        }

        [[nodiscard]] IndexType Inverse (IndexType it) const noexcept
        {
            return it;
        }
    };

    struct Swap02
    {
        [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 operator() (Dim3 i) const noexcept
        {
            return {i.z, i.y, i.x};
        }

        [[nodiscard]] AMREX_GPU_HOST_DEVICE Dim3 Inverse (Dim3 i) const noexcept
        {
            return {i.z, i.y, i.x};
        }

        [[nodiscard]] IndexType operator() (IndexType it) const noexcept
        {
            return it;
        }

        [[nodiscard]] IndexType Inverse (IndexType it) const noexcept
        {
            return it;
        }
    };

    // public for cuda
    template <typename F>
    void post_forward_doit (F const& post_forward);

private:

#if defined(AMREX_USE_CUDA)
    using FFTPlan = cufftHandle;
    using FFTPlan2 = FFTPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          cuComplex, cuDoubleComplex>;
#elif defined(AMREX_USE_HIP)
    using FFTPlan = rocfft_plan;
    using FFTPlan2 = FFTPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          float2, double2>;
#elif defined(AMREX_USE_SYCL)
    using FFTPlan = oneapi::mkl::dft::descriptor<
        std::is_same_v<float,T> ? oneapi::mkl::dft::precision::SINGLE
                                : oneapi::mkl::dft::precision::DOUBLE,
        oneapi::mkl::dft::domain::REAL> *;
    using FFTPlan2 = oneapi::mkl::dft::descriptor<
        std::is_same_v<float,T> ? oneapi::mkl::dft::precision::SINGLE
                                : oneapi::mkl::dft::precision::DOUBLE,
        oneapi::mkl::dft::domain::COMPLEX> *;
    using FFTComplex = GpuComplex<T>;
#else
    using FFTPlan = std::conditional_t<std::is_same_v<float,T>,
                                       fftwf_plan, fftw_plan>;
    using FFTPlan2 = FFTPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          fftwf_complex, fftw_complex>;
#endif

    void forward_doit (MF const& inmf, Scaling scaling = Scaling::none);
    void backward_doit (MF& outmf, Scaling scaling = Scaling::none);

    static void exec_r2c (FFTPlan plan, MF& in, cMF& out);
    static void exec_c2r (FFTPlan plan, cMF& in, MF& out);
    template <Direction direction>
    static void exec_c2c (FFTPlan2 plan, cMF& inout);

    template <typename P>
    static void destroy_plan (P plan);
    static std::pair<FFTPlan2,FFTPlan2> make_c2c_plans (cMF& inout);

    Box m_real_domain;
    Box m_spectral_domain_x;
    Box m_spectral_domain_y;
    Box m_spectral_domain_z;

    // assuming it's double for now
    FFTPlan m_fft_fwd_x;
    FFTPlan m_fft_bwd_x;
    FFTPlan2 m_fft_fwd_y;
    FFTPlan2 m_fft_bwd_y;
    FFTPlan2 m_fft_fwd_z;
    FFTPlan2 m_fft_bwd_z;

    // Comm meta-data. In the forward phase, we start with (x,y,z),
    // transpose to (y,x,z) and then (z,x,y). In the backward phase, we
    // perform inverse transpose.
    Swap01 m_dtos_x2y{};
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_x2y;
    //
    Swap01 m_dtos_y2x{};
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2x;
    //
    Swap02 m_dtos_y2z{};
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2z;
    //
    Swap02 m_dtos_z2y{};
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_z2y;

    // Optionally we need to copy from m_cz to user provided cMultiFab. xxxxx todo

    MF  m_rx;
    cMF m_cx;
    cMF m_cy;
    cMF m_cz;
};

template <typename T>
R2C<T>::R2C (Box const& domain)
    : m_real_domain(domain),
      m_spectral_domain_x(IntVect(0), IntVect(AMREX_D_DECL(domain.length(0)/2,
                                                           domain.bigEnd(1),
                                                           domain.bigEnd(2)))),
      m_spectral_domain_y(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(1),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(2)))),
      m_spectral_domain_z(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(2),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(1))))
{
    static_assert(std::is_same_v<float,T> || std::is_same_v<double,T>);
    AMREX_ALWAYS_ASSERT(m_real_domain.smallEnd() == 0 && m_real_domain.cellCentered());

    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();

    // xxxxx todo: need to handle cases there are more processes than 2d cells
    // xxxxx todo: 1d & 2d

    auto bax = amrex::decompose(m_real_domain, nprocs, {AMREX_D_DECL(false,true,true)});
    Vector<int> pmx(bax.size());
    std::iota(pmx.begin(), pmx.end(), 0);
    DistributionMapping dmx(std::move(pmx));
    m_rx.define(bax, dmx, 1, 0);

    {
        BoxList bl = bax.boxList();
        for (auto & b : bl) {
            b.setBig(0, m_spectral_domain_x.bigEnd(0));
        }
        BoxArray cbax(std::move(bl));
        m_cx.define(cbax, dmx, 1, 0);
    }

    // plans for x-direction
    {
        Box const local_box = m_rx.boxArray()[myproc];
        int n = local_box.length(0);
        int howmany = local_box.length(1) * local_box.length(2);

#if defined(AMREX_USE_CUDA)

        cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_R2C : CUFFT_D2Z;
        cufftType bwd_type = std::is_same_v<float,T> ? CUFFT_C2R : CUFFT_Z2D;
        AMREX_CUFFT_SAFE_CALL
            (cufftPlanMany(&m_fft_fwd_x, 1, &n, nullptr, 1, m_real_domain.length(0),
                           nullptr, 1, m_spectral_domain_x.length(0),
                           fwd_type, howmany));
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(m_fft_fwd_x, Gpu::gpuStream()));
        AMREX_CUFFT_SAFE_CALL
            (cufftPlanMany(&m_fft_bwd_x, 1, &n, nullptr, 1, m_spectral_domain_x.length(0),
                           nullptr, 1, m_real_domain.length(0),
                           bwd_type, howmany));
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(m_fft_bwd_x, Gpu::gpuStream()));

#elif defined(AMREX_USE_HIP)

        auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
        const std::size_t length = n;
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&m_fft_fwd_x, rocfft_placement_notinplace,
                                rocfft_transform_type_real_forward, prec, 1, &length, howmany,
                                nullptr));
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&m_fft_bwd_x, rocfft_placement_notinplace,
                                rocfft_transform_type_real_inverse, prec, 1, &length, howmany,
                                nullptr));

#elif defined(AMREX_USE_SYCL)

        m_fft_fwd_x = new std::remove_pointer_t<FFTPlan>(n);
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                               DFTI_NOT_INPLACE);
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                               howmany);
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE,
                               m_real_domain.length(0));
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE,
                               m_spectral_domain_x.length(0));
        std::array<std::int64_t,2> strides{0,1};
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES,
                               strides.data());
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES,
                               strides.data());
        m_fft_fwd_x->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                               oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        m_fft_fwd_x->commit(amrex::Gpu::Device::streamQueue());

        m_fft_bwd_x = m_fft_fwd_x;

#else /* FFTW */

        auto* in = m_rx[myproc].dataPtr();
        auto* out = (FFTComplex*)(m_cx[myproc].dataPtr());

        if constexpr (std::is_same_v<float,T>) {
            m_fft_fwd_x = fftwf_plan_many_dft_r2c
                (1, &n, howmany, in, nullptr, 1, m_real_domain.length(0),
                 out, nullptr, 1, m_spectral_domain_x.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

            m_fft_bwd_x = fftwf_plan_many_dft_c2r
                (1, &n, howmany, out, nullptr, 1, m_spectral_domain_x.length(0),
                 in, nullptr, 1, m_real_domain.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        } else {
            m_fft_fwd_x = fftw_plan_many_dft_r2c
                (1, &n, howmany, in, nullptr, 1, m_real_domain.length(0),
                 out, nullptr, 1, m_spectral_domain_x.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

            m_fft_bwd_x = fftw_plan_many_dft_c2r
                (1, &n, howmany, out, nullptr, 1, m_spectral_domain_x.length(0),
                 in, nullptr, 1, m_real_domain.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        }
#endif
    }

    auto cbay = amrex::decompose(m_spectral_domain_y, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping const& cdmy = dmx; // xxxxx todo
    m_cy.define(cbay, cdmy, 1, 0);

    std::tie(m_fft_fwd_y, m_fft_bwd_y) = make_c2c_plans(m_cy);

    // comm meta-data between x and y phases
    m_cmd_x2y = std::make_unique<MultiBlockCommMetaData>
        (m_cy.boxArray(), m_cy.DistributionMap(), m_spectral_domain_y,
         m_cx.boxArray(), m_cx.DistributionMap(), IntVect(0), m_dtos_x2y);
    m_cmd_y2x = std::make_unique<MultiBlockCommMetaData>
        (m_cx.boxArray(), m_cx.DistributionMap(), m_spectral_domain_x,
         m_cy.boxArray(), m_cy.DistributionMap(), IntVect(0), m_dtos_y2x);

    auto cbaz = amrex::decompose(m_spectral_domain_z, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping const& cdmz = dmx; // xxxxx todo
    m_cz.define(cbaz, cdmz, 1, 0);

    std::tie(m_fft_fwd_z, m_fft_bwd_z) = make_c2c_plans(m_cz);

    // comm meta-data between y and z phases
    m_cmd_y2z = std::make_unique<MultiBlockCommMetaData>
        (m_cz.boxArray(), m_cz.DistributionMap(), m_spectral_domain_z,
         m_cy.boxArray(), m_cy.DistributionMap(), IntVect(0), m_dtos_y2z);
    m_cmd_z2y = std::make_unique<MultiBlockCommMetaData>
        (m_cy.boxArray(), m_cy.DistributionMap(), m_spectral_domain_y,
         m_cz.boxArray(), m_cz.DistributionMap(), IntVect(0), m_dtos_z2y);
}

template <typename T>
template <typename P>
void R2C<T>::destroy_plan (P plan)
{
#if defined(AMREX_USE_CUDA)
    AMREX_CUFFT_SAFE_CALL(cufftDestroy(plan));
#elif defined(AMREX_USE_HIP)
    AMREX_ROCFFT_SAFE_CALL(rocfft_plan_destroy(plan));
#elif defined(AMREX_USE_SYCL)
    delete plan;
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_destroy_plan(plan);
    } else {
        fftw_destroy_plan(plan);
    }
#endif
}

template <typename T>
R2C<T>::~R2C<T> ()
{
    destroy_plan(m_fft_fwd_x);
    destroy_plan(m_fft_fwd_y);
    destroy_plan(m_fft_fwd_z);
#if !defined(AMREX_USE_SYCL)
    // For SYCL the same plans are used for both forward and backward transforms.
    destroy_plan(m_fft_bwd_x);
    destroy_plan(m_fft_bwd_y);
    destroy_plan(m_fft_bwd_z);
#endif
}

#ifdef AMREX_USE_HIP
namespace detail { void hip_execute (rocfft_plan plan, void **in, void **out); }
#endif

#ifdef AMREX_USE_SYCL
namespace detail
{
template <typename T, Direction direction, typename P, typename TI, typename TO>
void sycl_execute (P plan, TI* in, TO* out)
{
    std::size_t workspaceSize = 0;
    plan->get_value(oneapi::mkl::dft::config_param::WORKSPACE_BYTES,
                    &workspaceSize);
    auto* buffer = (T*)amrex::The_Arena()->alloc(workspaceSize);
    plan->set_workspace(buffer);
    sycl::event r;
    if (std::is_same_v<TI,TO>) {
        amrex::ignore_unused(in);
        if constexpr (direction == Direction::forward) {
            r = oneapi::mkl::dft::compute_forward(*plan, out);
        } else {
            r = oneapi::mkl::dft::compute_backward(*plan, out);
        }
    } else {
        if constexpr (direction == Direction::forward) {
            r = oneapi::mkl::dft::compute_forward(*plan, in, out);
        } else {
            r = oneapi::mkl::dft::compute_backward(*plan, in, out);
        }
    }
    r.wait();
    amrex::The_Arena()->free(buffer);
}
}
#endif

template <typename T>
void R2C<T>::exec_r2c (FFTPlan plan, MF& in, cMF& out)
{
#if defined(AMREX_USE_GPU)
    auto* pin = in[ParallelDescriptor::MyProc()].dataPtr();
    auto* pout = out[ParallelDescriptor::MyProc()].dataPtr();
#else
    amrex::ignore_unused(in,out);
#endif

#if defined(AMREX_USE_CUDA)
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecR2C(plan, pin, (FFTComplex*)pout));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecD2Z(plan, pin, (FFTComplex*)pout));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan, (void**)&pin, (void**)&pout);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,Direction::forward>(plan, pin, (std::complex<T>*)pout);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan);
    } else {
        fftw_execute(plan);
    }
#endif
}

template <typename T>
void R2C<T>::exec_c2r (FFTPlan plan, cMF& in, MF& out)
{
#if defined(AMREX_USE_GPU)
    auto* pin = in[ParallelDescriptor::MyProc()].dataPtr();
    auto* pout = out[ParallelDescriptor::MyProc()].dataPtr();
#else
    amrex::ignore_unused(in,out);
#endif

#if defined(AMREX_USE_CUDA)
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecC2R(plan, (FFTComplex*)pin, pout));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecZ2D(plan, (FFTComplex*)pin, pout));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan, (void**)&pin, (void**)&pout);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,Direction::backward>(plan, (std::complex<T>*)pin, pout);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan);
    } else {
        fftw_execute(plan);
    }
#endif
}

template <typename T>
template <Direction direction>
void R2C<T>::exec_c2c (FFTPlan2 plan, cMF& inout)
{
    amrex::ignore_unused(inout);
#if defined(AMREX_USE_GPU)
    auto* p = inout[ParallelDescriptor::MyProc()].dataPtr();
#endif

#if defined(AMREX_USE_CUDA)
    auto cufft_direction = (direction == Direction::forward) ? CUFFT_FORWARD : CUFFT_INVERSE;
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecC2C(plan, (FFTComplex*)p, (FFTComplex*)p,
                                           cufft_direction));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecZ2Z(plan, (FFTComplex*)p, (FFTComplex*)p,
                                           cufft_direction));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan, (void**)&p, (void**)&p);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,direction>(plan, (std::complex<T>*)p, (std::complex<T>*)p);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan);
    } else {
        fftw_execute(plan);
    }
#endif
}

template <typename T>
void R2C<T>::forward_doit (MF const& inmf, Scaling /*scaling*/)
{
    m_rx.ParallelCopy(inmf, 0, 0, 1);
    exec_r2c(m_fft_fwd_x, m_rx, m_cx);

    ParallelCopy(m_cy, m_cx, *m_cmd_x2y, 0, 0, 1, m_dtos_x2y);
    exec_c2c<Direction::forward>(m_fft_fwd_y, m_cy);

    ParallelCopy(m_cz, m_cy, *m_cmd_y2z, 0, 0, 1, m_dtos_y2z);
    exec_c2c<Direction::forward>(m_fft_fwd_z, m_cz);
}

template <typename T>
void R2C<T>::backward_doit (MF& outmf, Scaling /*scaling*/)
{
    exec_c2c<Direction::backward>(m_fft_bwd_z, m_cz);
    ParallelCopy(m_cy, m_cz, *m_cmd_z2y, 0, 0, 1, m_dtos_z2y);

    exec_c2c<Direction::backward>(m_fft_bwd_y, m_cy);
    ParallelCopy(m_cx, m_cy, *m_cmd_y2x, 0, 0, 1, m_dtos_y2x);

    exec_c2r(m_fft_bwd_x, m_cx, m_rx);
    outmf.ParallelCopy(m_rx, 0, 0, 1);
}

template <typename T>
std::pair<typename R2C<T>::FFTPlan2, typename R2C<T>::FFTPlan2>
R2C<T>::make_c2c_plans (cMF& inout)
{
    auto& fab = inout[ParallelDescriptor::MyProc()];
    Box const& local_box = fab.box();

    int n = local_box.length(0);
    int howmany = local_box.length(1) * local_box.length(2);

    FFTPlan2 fwd;
    FFTPlan2 bwd;

#if defined(AMREX_USE_CUDA)

    cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_C2C : CUFFT_Z2Z;
    cufftType bwd_type = std::is_same_v<float,T> ? CUFFT_C2C : CUFFT_Z2Z;
    AMREX_CUFFT_SAFE_CALL
        (cufftPlanMany(&fwd, 1, &n, nullptr, 1, n, nullptr, 1, n, fwd_type, howmany));
    AMREX_CUFFT_SAFE_CALL(cufftSetStream(fwd, Gpu::gpuStream()));
    AMREX_CUFFT_SAFE_CALL
        (cufftPlanMany(&bwd, 1, &n, nullptr, 1, n, nullptr, 1, n, bwd_type, howmany));
    AMREX_CUFFT_SAFE_CALL(cufftSetStream(bwd, Gpu::gpuStream()));

#elif defined(AMREX_USE_HIP)

    auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
    const std::size_t length = n;
    AMREX_ROCFFT_SAFE_CALL
        (rocfft_plan_create(&fwd, rocfft_placement_inplace,
                            rocfft_transform_type_complex_forward, prec, 1, &length, howmany,
                            nullptr));
    AMREX_ROCFFT_SAFE_CALL
        (rocfft_plan_create(&bwd, rocfft_placement_inplace,
                            rocfft_transform_type_complex_inverse, prec, 1, &length, howmany,
                            nullptr));

#elif defined(AMREX_USE_SYCL)

        fwd = new std::remove_pointer_t<FFTPlan2>(n);
        fwd->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                       DFTI_INPLACE);
        fwd->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                       howmany);
        fwd->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, n);
        fwd->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, n);
        std::array<std::int64_t,2> strides{0,1};
        fwd->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
        fwd->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides.data());
        fwd->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                               oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        fwd->commit(amrex::Gpu::Device::streamQueue());

        bwd = fwd;

#else
    auto* pinout = (FFTComplex*)fab.dataPtr();

    if constexpr (std::is_same_v<float,T>) {
        fwd = fftwf_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                  pinout, nullptr, 1, n, -1, FFTW_ESTIMATE);
        bwd = fftwf_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                  pinout, nullptr, 1, n, +1, FFTW_ESTIMATE);
    } else {
        fwd = fftw_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                 pinout, nullptr, 1, n, -1, FFTW_ESTIMATE);
        bwd = fftw_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                 pinout, nullptr, 1, n, +1, FFTW_ESTIMATE);
    }
#endif

    return {fwd,bwd};
}

template <typename T>
template <typename F>
void R2C<T>::post_forward_doit (F const& post_forward)
{
    auto& spectral_fab = m_cz[ParallelDescriptor::MyProc()];
    auto const& a = spectral_fab.array(); // m_cz's ordering is z,x,y
    ParallelFor(spectral_fab.box(), [=] AMREX_GPU_DEVICE (int iz, int jx, int ky)
    {
        post_forward(jx,ky,iz,a(iz,jx,ky));
    });
}

}

#endif
