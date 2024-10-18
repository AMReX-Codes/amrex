#ifndef AMREX_FFT_H_
#define AMREX_FFT_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_FFT_Helper.H>
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

template <typename T = Real>
class R2C
{
public:
    using MF = std::conditional_t<std::is_same_v<T,Real>,
                                  MultiFab, FabArray<BaseFab<T> > >;
    using cMF = FabArray<BaseFab<GpuComplex<T> > >;

    explicit R2C (Box const& domain, Info const& info = Info{});

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
    using VendorPlan = cufftHandle;
    using VendorPlan2 = VendorPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          cuComplex, cuDoubleComplex>;
#elif defined(AMREX_USE_HIP)
    using VendorPlan = rocfft_plan;
    using VendorPlan2 = VendorPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          float2, double2>;
#elif defined(AMREX_USE_SYCL)
    using VendorPlan = oneapi::mkl::dft::descriptor<
        std::is_same_v<float,T> ? oneapi::mkl::dft::precision::SINGLE
                                : oneapi::mkl::dft::precision::DOUBLE,
        oneapi::mkl::dft::domain::REAL> *;
    using VendorPlan2 = oneapi::mkl::dft::descriptor<
        std::is_same_v<float,T> ? oneapi::mkl::dft::precision::SINGLE
                                : oneapi::mkl::dft::precision::DOUBLE,
        oneapi::mkl::dft::domain::COMPLEX> *;
    using FFTComplex = GpuComplex<T>;
#else
    using VendorPlan = std::conditional_t<std::is_same_v<float,T>,
                                          fftwf_plan, fftw_plan>;
    using VendorPlan2 = VendorPlan;
    using FFTComplex = std::conditional_t<std::is_same_v<float,T>,
                                          fftwf_complex, fftw_complex>;
#endif

    struct Plan {
        bool defined = false;
        VendorPlan plan = 0; // NOLINT
    };

    struct Plan2 {
        bool defined = false;
        VendorPlan2 plan = 0; // NOLINT
    };

    template <typename FA>
    static typename FA::FABType::value_type *
    get_fab (FA& fa) {
        auto myproc = ParallelDescriptor::MyProc();
        if (myproc < fa.size()) {
            return fa.fabPtr(myproc);
        } else {
            return nullptr;
        }
    }

    void forward_doit (MF const& inmf, Scaling scaling = Scaling::none);
    void backward_doit (MF& outmf, Scaling scaling = Scaling::none);

    static void exec_r2c (Plan plan, MF& in, cMF& out);
    static void exec_c2r (Plan plan, cMF& in, MF& out);
    template <Direction direction>
    static void exec_c2c (Plan2 plan, cMF& inout);

    template <typename P>
    static void destroy_plan (P plan);
    static std::pair<Plan2,Plan2> make_c2c_plans (cMF& inout);

    Plan m_fft_fwd_x{};
    Plan m_fft_bwd_x{};
    Plan2 m_fft_fwd_y{};
    Plan2 m_fft_bwd_y{};
    Plan2 m_fft_fwd_z{};
    Plan2 m_fft_bwd_z{};

    // Comm meta-data. In the forward phase, we start with (x,y,z),
    // transpose to (y,x,z) and then (z,x,y). In the backward phase, we
    // perform inverse transpose.
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_x2y;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2x;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2z;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_z2y;
    Swap01 m_dtos_x2y{};
    Swap01 m_dtos_y2x{};
    Swap02 m_dtos_y2z{};
    Swap02 m_dtos_z2y{};

    // Optionally we need to copy from m_cz to user provided cMultiFab. xxxxx todo

    MF  m_rx;
    cMF m_cx;
    cMF m_cy;
    cMF m_cz;

    Box m_real_domain;
    Box m_spectral_domain_x;
    Box m_spectral_domain_y;
    Box m_spectral_domain_z;

    Info m_info;
};

template <typename T>
R2C<T>::R2C (Box const& domain, Info const& info)
    : m_real_domain(domain),
      m_spectral_domain_x(IntVect(0), IntVect(AMREX_D_DECL(domain.length(0)/2,
                                                           domain.bigEnd(1),
                                                           domain.bigEnd(2)))),
#if (AMREX_SPACEDIM >= 2)
      m_spectral_domain_y(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(1),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(2)))),
#if (AMREX_SPACEDIM == 3)
      m_spectral_domain_z(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(2),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(1)))),
#endif
#endif
      m_info(info)
{
    static_assert(std::is_same_v<float,T> || std::is_same_v<double,T>);
    AMREX_ALWAYS_ASSERT(m_real_domain.smallEnd() == 0 && m_real_domain.cellCentered());
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT(m_real_domain.length(2) > 1 || false == m_info.batch_mode);
#else
    AMREX_ALWAYS_ASSERT(false == m_info.batch_mode);
#endif

    int myproc = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();

    auto bax = amrex::decompose(m_real_domain, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping dmx;
    {
        Vector<int> pm(bax.size());
        std::iota(pm.begin(), pm.end(), 0);
        dmx.define(std::move(pm));
    }
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
    if (myproc < m_rx.size())
    {
        Box const local_box = m_rx.boxArray()[myproc];
        int n = local_box.length(0);
        int howmany = AMREX_D_TERM(1, *local_box.length(1), *local_box.length(2));

#if defined(AMREX_USE_CUDA)

        cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_R2C : CUFFT_D2Z;
        cufftType bwd_type = std::is_same_v<float,T> ? CUFFT_C2R : CUFFT_Z2D;
        AMREX_CUFFT_SAFE_CALL
            (cufftPlanMany(&m_fft_fwd_x.plan, 1, &n,
                           nullptr, 1, m_real_domain.length(0),
                           nullptr, 1, m_spectral_domain_x.length(0),
                           fwd_type, howmany));
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(m_fft_fwd_x.plan, Gpu::gpuStream()));
        AMREX_CUFFT_SAFE_CALL
            (cufftPlanMany(&m_fft_bwd_x.plan, 1, &n,
                           nullptr, 1, m_spectral_domain_x.length(0),
                           nullptr, 1, m_real_domain.length(0),
                           bwd_type, howmany));
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(m_fft_bwd_x.plan, Gpu::gpuStream()));

#elif defined(AMREX_USE_HIP)

        auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
        const std::size_t length = n;
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&m_fft_fwd_x.plan, rocfft_placement_notinplace,
                                rocfft_transform_type_real_forward, prec, 1,
                                &length, howmany, nullptr));
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&m_fft_bwd_x.plan, rocfft_placement_notinplace,
                                rocfft_transform_type_real_inverse, prec, 1,
                                &length, howmany, nullptr));

#elif defined(AMREX_USE_SYCL)

        m_fft_fwd_x.plan = new std::remove_pointer_t<VendorPlan>(n);
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                                    DFTI_NOT_INPLACE);
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                                    howmany);
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE,
                                    m_real_domain.length(0));
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE,
                                    m_spectral_domain_x.length(0));
        std::array<std::int64_t,2> strides{0,1};
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES,
                                    strides.data());
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES,
                                    strides.data());
        m_fft_fwd_x.plan->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                                    oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        m_fft_fwd_x.plan->commit(amrex::Gpu::Device::streamQueue());

        m_fft_bwd_x.plan = m_fft_fwd_x.plan;

#else /* FFTW */

        auto* in = m_rx[myproc].dataPtr();
        auto* out = (FFTComplex*)(m_cx[myproc].dataPtr());

        if constexpr (std::is_same_v<float,T>) {
            m_fft_fwd_x.plan = fftwf_plan_many_dft_r2c
                (1, &n, howmany, in, nullptr, 1, m_real_domain.length(0),
                 out, nullptr, 1, m_spectral_domain_x.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

            m_fft_bwd_x.plan = fftwf_plan_many_dft_c2r
                (1, &n, howmany, out, nullptr, 1, m_spectral_domain_x.length(0),
                 in, nullptr, 1, m_real_domain.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        } else {
            m_fft_fwd_x.plan = fftw_plan_many_dft_r2c
                (1, &n, howmany, in, nullptr, 1, m_real_domain.length(0),
                 out, nullptr, 1, m_spectral_domain_x.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);

            m_fft_bwd_x.plan = fftw_plan_many_dft_c2r
                (1, &n, howmany, out, nullptr, 1, m_spectral_domain_x.length(0),
                 in, nullptr, 1, m_real_domain.length(0),
                 FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
        }
#endif
        m_fft_fwd_x.defined = true;
        m_fft_bwd_x.defined = true;
    }

#if (AMREX_SPACEDIM >= 2)
    auto cbay = amrex::decompose(m_spectral_domain_y, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping cdmy;
    if (cbay.size() == dmx.size()) {
        cdmy = dmx;
    } else {
        Vector<int> pm(cbay.size());
        std::iota(pm.begin(), pm.end(), 0);
        cdmy.define(std::move(pm));
    }
    m_cy.define(cbay, cdmy, 1, 0);

    std::tie(m_fft_fwd_y, m_fft_bwd_y) = make_c2c_plans(m_cy);

    // comm meta-data between x and y phases
    m_cmd_x2y = std::make_unique<MultiBlockCommMetaData>
        (m_cy.boxArray(), m_cy.DistributionMap(), m_spectral_domain_y,
         m_cx.boxArray(), m_cx.DistributionMap(), IntVect(0), m_dtos_x2y);
    m_cmd_y2x = std::make_unique<MultiBlockCommMetaData>
        (m_cx.boxArray(), m_cx.DistributionMap(), m_spectral_domain_x,
         m_cy.boxArray(), m_cy.DistributionMap(), IntVect(0), m_dtos_y2x);

#if (AMREX_SPACEDIM == 3)
    if (false == m_info.batch_mode || m_real_domain.length(2) > 1) {
        auto cbaz = amrex::decompose(m_spectral_domain_z, nprocs, {false,true,true});
        DistributionMapping cdmz;
        if (cbaz.size() == dmx.size()) {
            cdmz = dmx;
        } else if (cbaz.size() == cdmy.size()) {
            cdmz = cdmy;
        } else {
            Vector<int> pm(cbaz.size());
            std::iota(pm.begin(), pm.end(), 0);
            cdmz.define(std::move(pm));
        }
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
#endif
#endif
}

template <typename T>
template <typename P>
void R2C<T>::destroy_plan (P plan)
{
    if (! plan.defined) { return; }

#if defined(AMREX_USE_CUDA)
    AMREX_CUFFT_SAFE_CALL(cufftDestroy(plan.plan));
#elif defined(AMREX_USE_HIP)
    AMREX_ROCFFT_SAFE_CALL(rocfft_plan_destroy(plan.plan));
#elif defined(AMREX_USE_SYCL)
    delete plan.plan;
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_destroy_plan(plan.plan);
    } else {
        fftw_destroy_plan(plan.plan);
    }
#endif

    plan.defined = false;
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
void R2C<T>::exec_r2c (Plan plan, MF& in, cMF& out)
{
    if (! plan.defined) { return; }

#if defined(AMREX_USE_GPU)
    auto* pin = in[ParallelDescriptor::MyProc()].dataPtr();
    auto* pout = out[ParallelDescriptor::MyProc()].dataPtr();
#else
    amrex::ignore_unused(in,out);
#endif

#if defined(AMREX_USE_CUDA)
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecR2C(plan.plan, pin, (FFTComplex*)pout));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecD2Z(plan.plan, pin, (FFTComplex*)pout));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan.plan, (void**)&pin, (void**)&pout);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,Direction::forward>(plan.plan, pin, (std::complex<T>*)pout);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan.plan);
    } else {
        fftw_execute(plan.plan);
    }
#endif
}

template <typename T>
void R2C<T>::exec_c2r (Plan plan, cMF& in, MF& out)
{
    if (! plan.defined) { return; }

#if defined(AMREX_USE_GPU)
    auto* pin = in[ParallelDescriptor::MyProc()].dataPtr();
    auto* pout = out[ParallelDescriptor::MyProc()].dataPtr();
#else
    amrex::ignore_unused(in,out);
#endif

#if defined(AMREX_USE_CUDA)
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecC2R(plan.plan, (FFTComplex*)pin, pout));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecZ2D(plan.plan, (FFTComplex*)pin, pout));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan.plan, (void**)&pin, (void**)&pout);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,Direction::backward>(plan.plan, (std::complex<T>*)pin, pout);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan.plan);
    } else {
        fftw_execute(plan.plan);
    }
#endif
}

template <typename T>
template <Direction direction>
void R2C<T>::exec_c2c (Plan2 plan, cMF& inout)
{
    if (! plan.defined) { return; }

    amrex::ignore_unused(inout);
#if defined(AMREX_USE_GPU)
    auto* p = inout[ParallelDescriptor::MyProc()].dataPtr();
#endif

#if defined(AMREX_USE_CUDA)
    auto cufft_direction = (direction == Direction::forward) ? CUFFT_FORWARD : CUFFT_INVERSE;
    if constexpr (std::is_same_v<float,T>) {
        AMREX_CUFFT_SAFE_CALL(cufftExecC2C(plan.plan, (FFTComplex*)p, (FFTComplex*)p,
                                           cufft_direction));
    } else {
        AMREX_CUFFT_SAFE_CALL(cufftExecZ2Z(plan.plan, (FFTComplex*)p, (FFTComplex*)p,
                                           cufft_direction));
    }
#elif defined(AMREX_USE_HIP)
    detail::hip_execute(plan.plan, (void**)&p, (void**)&p);
#elif defined(AMREX_USE_SYCL)
    detail::sycl_execute<T,direction>(plan.plan, (std::complex<T>*)p, (std::complex<T>*)p);
#else
    if constexpr (std::is_same_v<float,T>) {
        fftwf_execute(plan.plan);
    } else {
        fftw_execute(plan.plan);
    }
#endif
}

template <typename T>
void R2C<T>::forward_doit (MF const& inmf, Scaling /*scaling*/)
{
    m_rx.ParallelCopy(inmf, 0, 0, 1);
    exec_r2c(m_fft_fwd_x, m_rx, m_cx);

    if (                          m_cmd_x2y) {
        ParallelCopy(m_cy, m_cx, *m_cmd_x2y, 0, 0, 1, m_dtos_x2y);
    }
    exec_c2c<Direction::forward>(m_fft_fwd_y, m_cy);

    if (                          m_cmd_y2z) {
        ParallelCopy(m_cz, m_cy, *m_cmd_y2z, 0, 0, 1, m_dtos_y2z);
    }
    exec_c2c<Direction::forward>(m_fft_fwd_z, m_cz);
}

template <typename T>
void R2C<T>::backward_doit (MF& outmf, Scaling /*scaling*/)
{
    exec_c2c<Direction::backward>(m_fft_bwd_z, m_cz);
    if (                          m_cmd_z2y) {
        ParallelCopy(m_cy, m_cz, *m_cmd_z2y, 0, 0, 1, m_dtos_z2y);
    }

    exec_c2c<Direction::backward>(m_fft_bwd_y, m_cy);
    if (                          m_cmd_y2x) {
        ParallelCopy(m_cx, m_cy, *m_cmd_y2x, 0, 0, 1, m_dtos_y2x);
    }

    exec_c2r(m_fft_bwd_x, m_cx, m_rx);
    outmf.ParallelCopy(m_rx, 0, 0, 1);
}

template <typename T>
std::pair<typename R2C<T>::Plan2, typename R2C<T>::Plan2>
R2C<T>::make_c2c_plans (cMF& inout)
{
    Plan2 fwd;
    Plan2 bwd;

    auto* fab = get_fab(inout);
    if (!fab) { return {fwd, bwd};}

    Box const& local_box = fab->box();

    int n = local_box.length(0);
    int howmany = AMREX_D_TERM(1, *local_box.length(1), *local_box.length(2));

#if defined(AMREX_USE_CUDA)

    cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_C2C : CUFFT_Z2Z;
    cufftType bwd_type = std::is_same_v<float,T> ? CUFFT_C2C : CUFFT_Z2Z;
    AMREX_CUFFT_SAFE_CALL
        (cufftPlanMany(&fwd.plan, 1, &n, nullptr, 1, n, nullptr, 1, n, fwd_type, howmany));
    AMREX_CUFFT_SAFE_CALL(cufftSetStream(fwd.plan, Gpu::gpuStream()));
    AMREX_CUFFT_SAFE_CALL
        (cufftPlanMany(&bwd.plan, 1, &n, nullptr, 1, n, nullptr, 1, n, bwd_type, howmany));
    AMREX_CUFFT_SAFE_CALL(cufftSetStream(bwd.plan, Gpu::gpuStream()));

#elif defined(AMREX_USE_HIP)

    auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
    const std::size_t length = n;
    AMREX_ROCFFT_SAFE_CALL
        (rocfft_plan_create(&fwd.plan, rocfft_placement_inplace,
                            rocfft_transform_type_complex_forward, prec, 1, &length, howmany,
                            nullptr));
    AMREX_ROCFFT_SAFE_CALL
        (rocfft_plan_create(&bwd.plan, rocfft_placement_inplace,
                            rocfft_transform_type_complex_inverse, prec, 1, &length, howmany,
                            nullptr));

#elif defined(AMREX_USE_SYCL)

    fwd.plan = new std::remove_pointer_t<VendorPlan2>(n);
    fwd.plan->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                        DFTI_INPLACE);
    fwd.plan->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS,
                        howmany);
    fwd.plan->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, n);
    fwd.plan->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, n);
    std::array<std::int64_t,2> strides{0,1};
    fwd.plan->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
    fwd.plan->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides.data());
    fwd.plan->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                        oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
    fwd.plan->commit(amrex::Gpu::Device::streamQueue());

    bwd.plan = fwd.plan;

#else
    auto* pinout = (FFTComplex*)fab->dataPtr();

    if constexpr (std::is_same_v<float,T>) {
        fwd.plan = fftwf_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                       pinout, nullptr, 1, n, -1, FFTW_ESTIMATE);
        bwd.plan = fftwf_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                       pinout, nullptr, 1, n, +1, FFTW_ESTIMATE);
    } else {
        fwd.plan = fftw_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                      pinout, nullptr, 1, n, -1, FFTW_ESTIMATE);
        bwd.plan = fftw_plan_many_dft(1, &n, howmany, pinout, nullptr, 1, n,
                                      pinout, nullptr, 1, n, +1, FFTW_ESTIMATE);
    }
#endif

    fwd.defined = true;
    bwd.defined = true;

    return {fwd,bwd};
}

template <typename T>
template <typename F>
void R2C<T>::post_forward_doit (F const& post_forward)
{
    if (m_info.batch_mode) {
        amrex::Abort("xxxxx todo: post_forward");
    } else {
        if (                           ! m_cz.empty()) {
            auto* spectral_fab = get_fab(m_cz);
            if (spectral_fab) {
                auto const& a = spectral_fab->array(); // m_cz's ordering is z,x,y
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int iz, int jx, int ky)
                {
                    post_forward(jx,ky,iz,a(iz,jx,ky));
                });
            }
        } else if (                    ! m_cy.empty()) {
            auto* spectral_fab = get_fab(m_cy);
            if (spectral_fab) {
                auto const& a = spectral_fab->array(); // m_cy's ordering is y,x,z
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int iy, int jx, int k)
                {
                    post_forward(jx,iy,k,a(iy,jx,k));
                });
            }
        } else {
            auto* spectral_fab = get_fab(m_cx);
            if (spectral_fab) {
                auto const& a = spectral_fab->array();
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    post_forward(i,j,k,a(i,j,k));
                });
            }
        }
    }
}

}

#endif
