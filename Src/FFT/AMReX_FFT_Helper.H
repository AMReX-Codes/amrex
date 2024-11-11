#ifndef AMREX_FFT_HELPER_H_
#define AMREX_FFT_HELPER_H_
#include <AMReX_Config.H>

#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_DataAllocator.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_Enum.H>
#include <AMReX_Gpu.H>
#include <AMReX_GpuComplex.H>
#include <AMReX_Math.H>

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
#  if __has_include(<oneapi/mkl/dft.hpp>) // oneAPI 2025.0
#    include <oneapi/mkl/dft.hpp>
#else
#    define AMREX_USE_MKL_DFTI_2024 1
#    include <oneapi/mkl/dfti.hpp>
#  endif
#else
#  include <fftw3.h>
#endif

#include <algorithm>
#include <complex>
#include <limits>
#include <memory>
#include <tuple>
#include <utility>
#include <variant>

namespace amrex::FFT
{

enum struct Direction { forward, backward, both, none };

enum struct DomainStrategy { slab, pencil };

AMREX_ENUM( Boundary, periodic, even, odd );

enum struct Kind { none, r2c_f, r2c_b, c2c_f, c2c_b, r2r_ee_f, r2r_ee_b,
                   r2r_oo_f, r2r_oo_b, r2r_eo, r2r_oe };

struct Info
{
    //! Supported only in 3D. When batch_mode is true, FFT is performed on
    //! the first two dimensions only and the third dimension size is the
    //! batch size.
    bool batch_mode = false;

    //! Max number of processes to use
    int nprocs = std::numeric_limits<int>::max();

    Info& setBatchMode (bool x) { batch_mode = x; return *this; }
    Info& setNumProcs (int n) { nprocs = n; return *this; }
};

#ifdef AMREX_USE_HIP
namespace detail { void hip_execute (rocfft_plan plan, void **in, void **out); }
#endif

#ifdef AMREX_USE_SYCL
namespace detail
{
template <typename T, Direction direction, typename P, typename TI, typename TO>
void sycl_execute (P* plan, TI* in, TO* out)
{
#ifndef AMREX_USE_MKL_DFTI_2024
    std::int64_t workspaceSize = 0;
#else
    std::size_t workspaceSize = 0;
#endif
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
struct Plan
{
#if defined(AMREX_USE_CUDA)
    using VendorPlan = cufftHandle;
    using VendorComplex = std::conditional_t<std::is_same_v<float,T>,
                                             cuComplex, cuDoubleComplex>;
#elif defined(AMREX_USE_HIP)
    using VendorPlan = rocfft_plan;
    using VendorComplex = std::conditional_t<std::is_same_v<float,T>,
                                             float2, double2>;
#elif defined(AMREX_USE_SYCL)
    using mkl_desc_r = oneapi::mkl::dft::descriptor<std::is_same_v<float,T>
                                     ? oneapi::mkl::dft::precision::SINGLE
                                     : oneapi::mkl::dft::precision::DOUBLE,
                                     oneapi::mkl::dft::domain::REAL>;
    using mkl_desc_c = oneapi::mkl::dft::descriptor<std::is_same_v<float,T>
                                     ? oneapi::mkl::dft::precision::SINGLE
                                     : oneapi::mkl::dft::precision::DOUBLE,
                                     oneapi::mkl::dft::domain::COMPLEX>;
    using VendorPlan = std::variant<mkl_desc_r*,mkl_desc_c*>;
    using VendorComplex = std::complex<T>;
#else
    using VendorPlan = std::conditional_t<std::is_same_v<float,T>,
                                          fftwf_plan, fftw_plan>;
    using VendorComplex = std::conditional_t<std::is_same_v<float,T>,
                                             fftwf_complex, fftw_complex>;
#endif

    int n = 0;
    int howmany = 0;
    Kind kind = Kind::none;
    bool r2r_data_is_complex = false;
    bool defined = false;
    bool defined2 = false;
    VendorPlan plan{};
    VendorPlan plan2{};
    void* pf = nullptr;
    void* pb = nullptr;

#ifdef AMREX_USE_GPU
    void set_ptrs (void* p0, void* p1) {
        pf = p0;
        pb = p1;
    }
#endif

    void destroy ()
    {
        if (defined) {
            destroy_vendor_plan(plan);
            defined = false;
        }
#if !defined(AMREX_USE_GPU)
        if (defined2) {
            destroy_vendor_plan(plan2);
            defined2 = false;
        }
#endif
    }

    template <Direction D>
    void init_r2c (Box const& box, T* pr, VendorComplex* pc, bool is_2d_transform = false)
    {
        static_assert(D == Direction::forward || D == Direction::backward);

        int rank = is_2d_transform ? 2 : 1;

        kind = (D == Direction::forward) ? Kind::r2c_f : Kind::r2c_b;
        defined = true;
        pf = (void*)pr;
        pb = (void*)pc;

        int len[2] = {};
        if (rank == 1) {
            len[0] = box.length(0);
            len[1] = box.length(0); // Not used except for HIP. Yes it's `(0)`.
        } else {
            len[0] = box.length(1); // Most FFT libraries assume row-major ordering
            len[1] = box.length(0); // except for rocfft
        }
        int nr = (rank == 1) ? len[0] : len[0]*len[1];
        n = nr;
        int nc = (rank == 1) ? (len[0]/2+1) : (len[1]/2+1)*len[0];
#if (AMREX_SPACEDIM == 1)
        howmany = 1;
#else
        howmany = (rank == 1) ? AMREX_D_TERM(1, *box.length(1), *box.length(2))
                              : AMREX_D_TERM(1, *1            , *box.length(2));
#endif

        amrex::ignore_unused(nc);

#if defined(AMREX_USE_CUDA)

        AMREX_CUFFT_SAFE_CALL(cufftCreate(&plan));
        AMREX_CUFFT_SAFE_CALL(cufftSetAutoAllocation(plan, 0));
        std::size_t work_size;
        if constexpr (D == Direction::forward) {
            cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_R2C : CUFFT_D2Z;
            AMREX_CUFFT_SAFE_CALL
                (cufftMakePlanMany(plan, rank, len, nullptr, 1, nr, nullptr, 1, nc, fwd_type, howmany, &work_size));
        } else {
            cufftType bwd_type = std::is_same_v<float,T> ? CUFFT_C2R : CUFFT_Z2D;
            AMREX_CUFFT_SAFE_CALL
                (cufftMakePlanMany(plan, rank, len, nullptr, 1, nc, nullptr, 1, nr, bwd_type, howmany, &work_size));
        }

#elif defined(AMREX_USE_HIP)

        auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
        // switch to column-major ordering
        std::size_t length[2] = {std::size_t(len[1]), std::size_t(len[0])};
        if constexpr (D == Direction::forward) {
            AMREX_ROCFFT_SAFE_CALL
                (rocfft_plan_create(&plan, rocfft_placement_notinplace,
                                    rocfft_transform_type_real_forward, prec, rank,
                                    length, howmany, nullptr));
        } else {
            AMREX_ROCFFT_SAFE_CALL
                (rocfft_plan_create(&plan, rocfft_placement_notinplace,
                                    rocfft_transform_type_real_inverse, prec, rank,
                                    length, howmany, nullptr));
        }

#elif defined(AMREX_USE_SYCL)

        mkl_desc_r* pp;
        if (rank == 1) {
            pp = new mkl_desc_r(len[0]);
        } else {
            pp = new mkl_desc_r({std::int64_t(len[0]), std::int64_t(len[1])});
        }
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                      oneapi::mkl::dft::config_value::NOT_INPLACE);
#else
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_NOT_INPLACE);
#endif
        pp->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, howmany);
        pp->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nr);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nc);
        std::vector<std::int64_t> strides;
        strides.push_back(0);
        if (rank == 2) { strides.push_back(len[1]); }
        strides.push_back(1);
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides);
#else
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides.data());
#endif
        pp->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                      oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        pp->commit(amrex::Gpu::Device::streamQueue());
        plan = pp;

#else /* FFTW */

        if constexpr (std::is_same_v<float,T>) {
            if constexpr (D == Direction::forward) {
                plan = fftwf_plan_many_dft_r2c
                    (rank, len, howmany, pr, nullptr, 1, nr, pc, nullptr, 1, nc,
                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            } else {
                plan = fftwf_plan_many_dft_c2r
                    (rank, len, howmany, pc, nullptr, 1, nc, pr, nullptr, 1, nr,
                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }
        } else {
            if constexpr (D == Direction::forward) {
                plan = fftw_plan_many_dft_r2c
                    (rank, len, howmany, pr, nullptr, 1, nr, pc, nullptr, 1, nc,
                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            } else {
                plan = fftw_plan_many_dft_c2r
                    (rank, len, howmany, pc, nullptr, 1, nc, pr, nullptr, 1, nr,
                     FFTW_ESTIMATE | FFTW_DESTROY_INPUT);
            }
        }
#endif
    }

    template <Direction D, int M>
    void init_r2c (IntVectND<M> const& fft_size, void*, void*, bool cache);

    template <Direction D>
    void init_c2c (Box const& box, VendorComplex* p)
    {
        static_assert(D == Direction::forward || D == Direction::backward);

        kind = (D == Direction::forward) ? Kind::c2c_f : Kind::c2c_b;
        defined = true;
        pf = (void*)p;
        pb = (void*)p;

        n = box.length(0);
        howmany = AMREX_D_TERM(1, *box.length(1), *box.length(2));

#if defined(AMREX_USE_CUDA)
        AMREX_CUFFT_SAFE_CALL(cufftCreate(&plan));
        AMREX_CUFFT_SAFE_CALL(cufftSetAutoAllocation(plan, 0));

        cufftType t = std::is_same_v<float,T> ? CUFFT_C2C : CUFFT_Z2Z;
        std::size_t work_size;
        AMREX_CUFFT_SAFE_CALL
            (cufftMakePlanMany(plan, 1, &n, nullptr, 1, n, nullptr, 1, n, t, howmany, &work_size));

#elif defined(AMREX_USE_HIP)

        auto prec = std::is_same_v<float,T> ? rocfft_precision_single
                                            : rocfft_precision_double;
        auto dir= (D == Direction::forward) ? rocfft_transform_type_complex_forward
                                            : rocfft_transform_type_complex_inverse;
        const std::size_t length = n;
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&plan, rocfft_placement_inplace, dir, prec, 1,
                                &length, howmany, nullptr));

#elif defined(AMREX_USE_SYCL)

        auto* pp = new mkl_desc_c(n);
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                      oneapi::mkl::dft::config_value::INPLACE);
#else
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_INPLACE);
#endif
        pp->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, howmany);
        pp->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, n);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, n);
        std::vector<std::int64_t> strides = {0,1};
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides);
#else
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides.data());
#endif
        pp->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                      oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        pp->commit(amrex::Gpu::Device::streamQueue());
        plan = pp;

#else /* FFTW */

        if constexpr (std::is_same_v<float,T>) {
            if constexpr (D == Direction::forward) {
                plan = fftwf_plan_many_dft
                    (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, -1,
                     FFTW_ESTIMATE);
            } else {
                plan = fftwf_plan_many_dft
                    (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, +1,
                     FFTW_ESTIMATE);
            }
        } else {
            if constexpr (D == Direction::forward) {
                plan = fftw_plan_many_dft
                    (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, -1,
                     FFTW_ESTIMATE);
            } else {
                plan = fftw_plan_many_dft
                    (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, +1,
                     FFTW_ESTIMATE);
            }
        }
#endif
    }

#ifndef AMREX_USE_GPU
    template <Direction D>
    fftw_r2r_kind get_fftw_kind (std::pair<Boundary,Boundary> const& bc)
    {
        if (bc.first == Boundary::even && bc.second == Boundary::even)
        {
            return (D == Direction::forward) ? FFTW_REDFT10 : FFTW_REDFT01;
        }
        else if (bc.first == Boundary::even && bc.second == Boundary::odd)
        {
            return FFTW_REDFT11;
        }
        else if (bc.first == Boundary::odd && bc.second == Boundary::even)
        {
            return FFTW_RODFT11;
        }
        else if (bc.first == Boundary::odd && bc.second == Boundary::odd)
        {
            return (D == Direction::forward) ? FFTW_RODFT10 : FFTW_RODFT01;
        }
        else {
            amrex::Abort("FFT: unsupported BC");
            return fftw_r2r_kind{};
        }

    }
#endif

    template <Direction D>
    Kind get_r2r_kind (std::pair<Boundary,Boundary> const& bc)
    {
        if (bc.first == Boundary::even && bc.second == Boundary::even)
        {
            return (D == Direction::forward) ? Kind::r2r_ee_f : Kind::r2r_ee_b;
        }
        else if (bc.first == Boundary::even && bc.second == Boundary::odd)
        {
            return Kind::r2r_eo;
        }
        else if (bc.first == Boundary::odd && bc.second == Boundary::even)
        {
            return Kind::r2r_oe;
        }
        else if (bc.first == Boundary::odd && bc.second == Boundary::odd)
        {
            return (D == Direction::forward) ? Kind::r2r_oo_f : Kind::r2r_oo_b;
        }
        else {
            amrex::Abort("FFT: unsupported BC");
            return Kind::none;
        }

    }

    template <Direction D>
    void init_r2r (Box const& box, T* p, std::pair<Boundary,Boundary> const& bc,
                   int howmany_initval = 1)
    {
        static_assert(D == Direction::forward || D == Direction::backward);

        kind = get_r2r_kind<D>(bc);
        defined = true;
        pf = (void*)p;
        pb = (void*)p;

        n = box.length(0);
        howmany = AMREX_D_TERM(howmany_initval, *box.length(1), *box.length(2));

#if defined(AMREX_USE_GPU)
        int nex=0;
        if (bc.first == Boundary::odd && bc.second == Boundary::odd &&
            Direction::forward == D) {
            nex = 2*n;
        } else if (bc.first == Boundary::odd && bc.second == Boundary::odd &&
            Direction::backward == D) {
            nex = 4*n;
        } else if (bc.first == Boundary::even && bc.second == Boundary::even &&
            Direction::forward == D) {
            nex = 2*n;
        } else if (bc.first == Boundary::even && bc.second == Boundary::even &&
            Direction::backward == D) {
            nex = 4*n;
        } else if ((bc.first == Boundary::even && bc.second == Boundary::odd) ||
                   (bc.first == Boundary::odd && bc.second == Boundary::even)) {
            nex = 4*n;
        } else {
            amrex::Abort("FFT: unsupported BC");
        }
        int nc = (nex/2) + 1;

#if defined (AMREX_USE_CUDA)

        AMREX_CUFFT_SAFE_CALL(cufftCreate(&plan));
        AMREX_CUFFT_SAFE_CALL(cufftSetAutoAllocation(plan, 0));
        cufftType fwd_type = std::is_same_v<float,T> ? CUFFT_R2C : CUFFT_D2Z;
        std::size_t work_size;
        AMREX_CUFFT_SAFE_CALL
            (cufftMakePlanMany(plan, 1, &nex, nullptr, 1, nc*2, nullptr, 1, nc, fwd_type, howmany, &work_size));

#elif defined(AMREX_USE_HIP)

        amrex::ignore_unused(nc);
        auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
        const std::size_t length = nex;
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&plan, rocfft_placement_inplace,
                                rocfft_transform_type_real_forward, prec, 1,
                                &length, howmany, nullptr));

#elif defined(AMREX_USE_SYCL)

        auto* pp = new mkl_desc_r(nex);
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                      oneapi::mkl::dft::config_value::INPLACE);
#else
        pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_INPLACE);
#endif
        pp->set_value(oneapi::mkl::dft::config_param::NUMBER_OF_TRANSFORMS, howmany);
        pp->set_value(oneapi::mkl::dft::config_param::FWD_DISTANCE, nc*2);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_DISTANCE, nc);
        std::vector<std::int64_t> strides = {0,1};
#ifndef AMREX_USE_MKL_DFTI_2024
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides);
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides);
#else
        pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
        pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides.data());
#endif
        pp->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                      oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
        pp->commit(amrex::Gpu::Device::streamQueue());
        plan = pp;

#endif

#else /* FFTW */
        auto fftw_kind = get_fftw_kind<D>(bc);
        if constexpr (std::is_same_v<float,T>) {
            plan = fftwf_plan_many_r2r
                (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, &fftw_kind,
                 FFTW_ESTIMATE);
        } else {
            plan = fftw_plan_many_r2r
                (1, &n, howmany, p, nullptr, 1, n, p, nullptr, 1, n, &fftw_kind,
                 FFTW_ESTIMATE);
        }
#endif
    }

    template <Direction D>
    void init_r2r (Box const& box, VendorComplex* pc,
                   std::pair<Boundary,Boundary> const& bc)
    {
        static_assert(D == Direction::forward || D == Direction::backward);

        auto* p = (T*)pc;

#if defined(AMREX_USE_GPU)

        init_r2r<D>(box, p, bc, 2);
        r2r_data_is_complex = true;

#else

        kind = get_r2r_kind<D>(bc);
        defined = true;
        pf = (void*)p;
        pb = (void*)p;

        n = box.length(0);
        howmany = AMREX_D_TERM(1, *box.length(1), *box.length(2));

        defined2 = true;
        auto fftw_kind = get_fftw_kind<D>(bc);
        if constexpr (std::is_same_v<float,T>) {
            plan = fftwf_plan_many_r2r
                (1, &n, howmany, p, nullptr, 2, n*2, p, nullptr, 2, n*2, &fftw_kind,
                 FFTW_ESTIMATE);
            plan2 = fftwf_plan_many_r2r
                (1, &n, howmany, p+1, nullptr, 2, n*2, p+1, nullptr, 2, n*2, &fftw_kind,
                 FFTW_ESTIMATE);
        } else {
            plan = fftw_plan_many_r2r
                (1, &n, howmany, p, nullptr, 2, n*2, p, nullptr, 2, n*2, &fftw_kind,
                 FFTW_ESTIMATE);
            plan2 = fftw_plan_many_r2r
                (1, &n, howmany, p+1, nullptr, 2, n*2, p+1, nullptr, 2, n*2, &fftw_kind,
                 FFTW_ESTIMATE);
        }
#endif
    }

    template <Direction D>
    void compute_r2c ()
    {
        static_assert(D == Direction::forward || D == Direction::backward);
        if (!defined) { return; }

        using TI = std::conditional_t<(D == Direction::forward), T, VendorComplex>;
        using TO = std::conditional_t<(D == Direction::backward), T, VendorComplex>;
        auto* pi = (TI*)((D == Direction::forward) ? pf : pb);
        auto* po = (TO*)((D == Direction::forward) ? pb : pf);

#if defined(AMREX_USE_CUDA)
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(plan, Gpu::gpuStream()));

        std::size_t work_size = 0;
        AMREX_CUFFT_SAFE_CALL(cufftGetSize(plan, &work_size));

        auto* work_area = The_Arena()->alloc(work_size);
        AMREX_CUFFT_SAFE_CALL(cufftSetWorkArea(plan, work_area));

        if constexpr (D == Direction::forward) {
            if constexpr (std::is_same_v<float,T>) {
                AMREX_CUFFT_SAFE_CALL(cufftExecR2C(plan, pi, po));
            } else {
                AMREX_CUFFT_SAFE_CALL(cufftExecD2Z(plan, pi, po));
            }
        } else {
            if constexpr (std::is_same_v<float,T>) {
                AMREX_CUFFT_SAFE_CALL(cufftExecC2R(plan, pi, po));
            } else {
                AMREX_CUFFT_SAFE_CALL(cufftExecZ2D(plan, pi, po));
            }
        }
        Gpu::streamSynchronize();
        The_Arena()->free(work_area);
#elif defined(AMREX_USE_HIP)
        detail::hip_execute(plan, (void**)&pi, (void**)&po);
#elif defined(AMREX_USE_SYCL)
        detail::sycl_execute<T,D>(std::get<0>(plan), pi, po);
#else
        amrex::ignore_unused(pi,po);
        if constexpr (std::is_same_v<float,T>) {
            fftwf_execute(plan);
        } else {
            fftw_execute(plan);
        }
#endif
    }

    template <Direction D>
    void compute_c2c ()
    {
        static_assert(D == Direction::forward || D == Direction::backward);
        if (!defined) { return; }

        auto* p = (VendorComplex*)pf;

#if defined(AMREX_USE_CUDA)
        AMREX_CUFFT_SAFE_CALL(cufftSetStream(plan, Gpu::gpuStream()));

        std::size_t work_size = 0;
        AMREX_CUFFT_SAFE_CALL(cufftGetSize(plan, &work_size));

        auto* work_area = The_Arena()->alloc(work_size);
        AMREX_CUFFT_SAFE_CALL(cufftSetWorkArea(plan, work_area));

        auto dir = (D == Direction::forward) ? CUFFT_FORWARD : CUFFT_INVERSE;
        if constexpr (std::is_same_v<float,T>) {
            AMREX_CUFFT_SAFE_CALL(cufftExecC2C(plan, p, p, dir));
        } else {
            AMREX_CUFFT_SAFE_CALL(cufftExecZ2Z(plan, p, p, dir));
        }
        Gpu::streamSynchronize();
        The_Arena()->free(work_area);
#elif defined(AMREX_USE_HIP)
        detail::hip_execute(plan, (void**)&p, (void**)&p);
#elif defined(AMREX_USE_SYCL)
        detail::sycl_execute<T,D>(std::get<1>(plan), p, p);
#else
        amrex::ignore_unused(p);
        if constexpr (std::is_same_v<float,T>) {
            fftwf_execute(plan);
        } else {
            fftw_execute(plan);
        }
#endif
    }

#ifdef AMREX_USE_GPU
    [[nodiscard]] void* alloc_scratch_space () const
    {
        int nc = 0;
        if (kind == Kind::r2r_oo_f || kind == Kind::r2r_ee_f) {
            nc = n + 1;
        } else if (kind == Kind::r2r_oo_b || kind == Kind::r2r_ee_b ||
                   kind == Kind::r2r_oe || kind == Kind::r2r_eo) {
            nc = 2*n+1;
        } else {
            amrex::Abort("FFT: alloc_scratch_space: unsupported kind");
        }
        return The_Arena()->alloc(sizeof(GpuComplex<T>)*nc*howmany);
    }

    static void free_scratch_space (void* p) { The_Arena()->free(p); }

    void pack_r2r_buffer (void* pbuf, T const* psrc) const
    {
        auto* pdst = (T*) pbuf;
        if (kind == Kind::r2r_oo_f || kind == Kind::r2r_ee_f) {
            T sign = (kind == Kind::r2r_oo_f) ? T(-1) : T(1);
            int ostride = (n+1)*2;
            int istride = n;
            int nex = 2*n;
            int norig = n;
            Long nelems = Long(nex)*howmany;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto* po = pdst + (2*batch+ir)*ostride + i;
                        auto const* pi = psrc + 2*batch*istride + ir;
                        if (i < norig) {
                            *po = pi[i*2];
                        } else {
                            *po = sign * pi[(2*norig-1-i)*2];
                        }
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    auto* po = pdst + batch*ostride + i;
                    auto const* pi = psrc + batch*istride;
                    if (i < norig) {
                        *po = pi[i];
                    } else {
                        *po = sign * pi[2*norig-1-i];
                    }
                });
            }
        } else if (kind == Kind::r2r_oo_b) {
            int ostride = (2*n+1)*2;
            int istride = n;
            int nex = 4*n;
            int norig = n;
            Long nelems = Long(nex)*howmany;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto* po = pdst + (2*batch+ir)*ostride + i;
                        auto const* pi = psrc + 2*batch*istride + ir;
                        if (i < norig) {
                            *po = pi[i*2];
                        } else if (i < (2*norig-1)) {
                            *po = pi[(2*norig-2-i)*2];
                        } else if (i == (2*norig-1)) {
                            *po = T(0);
                        } else if (i < (3*norig)) {
                            *po = -pi[(i-2*norig)*2];
                        } else if (i < (4*norig-1)) {
                            *po = -pi[(4*norig-2-i)*2];
                        } else {
                            *po = T(0);
                        }
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    auto* po = pdst + batch*ostride + i;
                    auto const* pi = psrc + batch*istride;
                    if (i < norig) {
                        *po = pi[i];
                    } else if (i < (2*norig-1)) {
                        *po = pi[2*norig-2-i];
                    } else if (i == (2*norig-1)) {
                        *po = T(0);
                    } else if (i < (3*norig)) {
                        *po = -pi[i-2*norig];
                    } else if (i < (4*norig-1)) {
                        *po = -pi[4*norig-2-i];
                    } else {
                        *po = T(0);
                    }
                });
            }
        } else if (kind == Kind::r2r_ee_b) {
            int ostride = (2*n+1)*2;
            int istride = n;
            int nex = 4*n;
            int norig = n;
            Long nelems = Long(nex)*howmany;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto* po = pdst + (2*batch+ir)*ostride + i;
                        auto const* pi = psrc + 2*batch*istride + ir;
                        if (i < norig) {
                            *po = pi[i*2];
                        } else if (i == norig) {
                            *po = T(0);
                        } else if (i < (2*norig+1)) {
                            *po = -pi[(2*norig-i)*2];
                        } else if (i < (3*norig)) {
                            *po = -pi[(i-2*norig)*2];
                        } else if (i == 3*norig) {
                            *po = T(0);
                        } else {
                            *po = pi[(4*norig-i)*2];
                        }
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    auto* po = pdst + batch*ostride + i;
                    auto const* pi = psrc + batch*istride;
                    if (i < norig) {
                        *po = pi[i];
                    } else if (i == norig) {
                        *po = T(0);
                    } else if (i < (2*norig+1)) {
                        *po = -pi[2*norig-i];
                    } else if (i < (3*norig)) {
                        *po = -pi[i-2*norig];
                    } else if (i == 3*norig) {
                        *po = T(0);
                    } else {
                        *po = pi[4*norig-i];
                    }
                });
            }
        } else if (kind == Kind::r2r_eo) {
            int ostride = (2*n+1)*2;
            int istride = n;
            int nex = 4*n;
            int norig = n;
            Long nelems = Long(nex)*howmany;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto* po = pdst + (2*batch+ir)*ostride + i;
                        auto const* pi = psrc + 2*batch*istride + ir;
                        if (i < norig) {
                            *po = pi[i*2];
                        } else if (i < (2*norig)) {
                            *po = -pi[(2*norig-1-i)*2];
                        } else if (i < (3*norig)) {
                            *po = -pi[(i-2*norig)*2];
                        } else {
                            *po = pi[(4*norig-1-i)*2];
                        }
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    auto* po = pdst + batch*ostride + i;
                    auto const* pi = psrc + batch*istride;
                    if (i < norig) {
                        *po = pi[i];
                    } else if (i < (2*norig)) {
                        *po = -pi[2*norig-1-i];
                    } else if (i < (3*norig)) {
                        *po = -pi[i-2*norig];
                    } else {
                        *po = pi[4*norig-1-i];
                    }
                });
            }
        } else if (kind == Kind::r2r_oe) {
            int ostride = (2*n+1)*2;
            int istride = n;
            int nex = 4*n;
            int norig = n;
            Long nelems = Long(nex)*howmany;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto* po = pdst + (2*batch+ir)*ostride + i;
                        auto const* pi = psrc + 2*batch*istride + ir;
                        if (i < norig) {
                            *po = pi[i*2];
                        } else if (i < (2*norig)) {
                            *po = pi[(2*norig-1-i)*2];
                        } else if (i < (3*norig)) {
                            *po = -pi[(i-2*norig)*2];
                        } else {
                            *po = -pi[(4*norig-1-i)*2];
                        }
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(nex);
                    auto i = int(ielem - batch*nex);
                    auto* po = pdst + batch*ostride + i;
                    auto const* pi = psrc + batch*istride;
                    if (i < norig) {
                        *po = pi[i];
                    } else if (i < (2*norig)) {
                        *po = pi[2*norig-1-i];
                    } else if (i < (3*norig)) {
                        *po = -pi[i-2*norig];
                    } else {
                        *po = -pi[4*norig-1-i];
                    }
                });
            }
        } else {
            amrex::Abort("FFT: pack_r2r_buffer: unsupported kind");
        }
    }

    void unpack_r2r_buffer (T* pdst, void const* pbuf) const
    {
        auto const* psrc = (GpuComplex<T> const*) pbuf;
        int norig = n;
        Long nelems = Long(norig)*howmany;
        int ostride = n;

        if (kind == Kind::r2r_oo_f) {
            int istride = n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(k+1)/T(2*norig));
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+k+1];
                        pdst[2*batch*ostride+ir+k*2] = s * yk.real() - c * yk.imag();
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(k+1)/T(2*norig));
                    auto const& yk = psrc[batch*istride+k+1];
                    pdst[batch*ostride+k] = s * yk.real() - c * yk.imag();
                });
            }
        } else if (kind == Kind::r2r_oo_b) {
            int istride = 2*n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(2*k+1)/T(2*norig));
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+2*k+1];
                        pdst[2*batch*ostride+ir+k*2] = T(0.5)*(s * yk.real() - c * yk.imag());
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(2*k+1)/T(2*norig));
                    auto const& yk = psrc[batch*istride+2*k+1];
                    pdst[batch*ostride+k] = T(0.5)*(s * yk.real() - c * yk.imag());
                });
            }
        } else if (kind == Kind::r2r_ee_f) {
            int istride = n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(k)/T(2*norig));
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+k];
                        pdst[2*batch*ostride+ir+k*2] = c * yk.real() + s * yk.imag();
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi(T(k)/T(2*norig));
                    auto const& yk = psrc[batch*istride+k];
                    pdst[batch*ostride+k] = c * yk.real() + s * yk.imag();
                });
            }
        } else if (kind == Kind::r2r_ee_b) {
            int istride = 2*n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+2*k+1];
                        pdst[2*batch*ostride+ir+k*2] = T(0.5) * yk.real();
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto const& yk = psrc[batch*istride+2*k+1];
                    pdst[batch*ostride+k] = T(0.5) * yk.real();
                });
            }
        } else if (kind == Kind::r2r_eo) {
            int istride = 2*n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi((k+T(0.5))/T(2*norig));
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+2*k+1];
                        pdst[2*batch*ostride+ir+k*2] = T(0.5) * (c * yk.real() + s * yk.imag());
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi((k+T(0.5))/T(2*norig));
                    auto const& yk = psrc[batch*istride+2*k+1];
                    pdst[batch*ostride+k] = T(0.5) * (c * yk.real() + s * yk.imag());
                });
            }
        } else if (kind == Kind::r2r_oe) {
            int istride = 2*n+1;
            if (r2r_data_is_complex) {
                ParallelFor(nelems/2, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi((k+T(0.5))/T(2*norig));
                    for (int ir = 0; ir < 2; ++ir) {
                        auto const& yk = psrc[(2*batch+ir)*istride+2*k+1];
                        pdst[2*batch*ostride+ir+k*2] = T(0.5) * (s * yk.real() - c * yk.imag());
                    }
                });
            } else {
                ParallelFor(nelems, [=] AMREX_GPU_DEVICE (Long ielem)
                {
                    auto batch = ielem / Long(norig);
                    auto k = int(ielem - batch*norig);
                    auto [s, c] = Math::sincospi((k+T(0.5))/T(2*norig));
                    auto const& yk = psrc[batch*istride+2*k+1];
                    pdst[batch*ostride+k] = T(0.5) * (s * yk.real() - c * yk.imag());
                });
            }
        } else {
            amrex::Abort("FFT: unpack_r2r_buffer: unsupported kind");
        }
    }
#endif

    template <Direction D>
    void compute_r2r ()
    {
        static_assert(D == Direction::forward || D == Direction::backward);
        if (!defined) { return; }

#if defined(AMREX_USE_GPU)

        auto* pscratch = alloc_scratch_space();

        pack_r2r_buffer(pscratch, (T*)((D == Direction::forward) ? pf : pb));

#if defined(AMREX_USE_CUDA)

        AMREX_CUFFT_SAFE_CALL(cufftSetStream(plan, Gpu::gpuStream()));

        std::size_t work_size = 0;
        AMREX_CUFFT_SAFE_CALL(cufftGetSize(plan, &work_size));

        auto* work_area = The_Arena()->alloc(work_size);
        AMREX_CUFFT_SAFE_CALL(cufftSetWorkArea(plan, work_area));

        if constexpr (std::is_same_v<float,T>) {
            AMREX_CUFFT_SAFE_CALL(cufftExecR2C(plan, (T*)pscratch, (VendorComplex*)pscratch));
        } else {
            AMREX_CUFFT_SAFE_CALL(cufftExecD2Z(plan, (T*)pscratch, (VendorComplex*)pscratch));
        }

#elif defined(AMREX_USE_HIP)
        detail::hip_execute(plan, (void**)&pscratch, (void**)&pscratch);
#elif defined(AMREX_USE_SYCL)
        detail::sycl_execute<T,Direction::forward>(std::get<0>(plan), (T*)pscratch, (VendorComplex*)pscratch);
#endif

        unpack_r2r_buffer((T*)((D == Direction::forward) ? pb : pf), pscratch);

        Gpu::streamSynchronize();
        free_scratch_space(pscratch);
#if defined(AMREX_USE_CUDA)
        The_Arena()->free(work_area);
#endif

#else /* FFTW */

        if constexpr (std::is_same_v<float,T>) {
            fftwf_execute(plan);
            if (defined2) { fftwf_execute(plan2); }
        } else {
            fftw_execute(plan);
            if (defined2) { fftw_execute(plan2); }
        }

#endif
    }

    static void destroy_vendor_plan (VendorPlan plan)
    {
#if defined(AMREX_USE_CUDA)
        AMREX_CUFFT_SAFE_CALL(cufftDestroy(plan));
#elif defined(AMREX_USE_HIP)
        AMREX_ROCFFT_SAFE_CALL(rocfft_plan_destroy(plan));
#elif defined(AMREX_USE_SYCL)
        std::visit([](auto&& p) { delete p; }, plan);
#else
        if constexpr (std::is_same_v<float,T>) {
            fftwf_destroy_plan(plan);
        } else {
            fftw_destroy_plan(plan);
        }
#endif
    }
};

using Key = std::tuple<IntVectND<3>,Direction,Kind>;
using PlanD = typename Plan<double>::VendorPlan;
using PlanF = typename Plan<float>::VendorPlan;

PlanD* get_vendor_plan_d (Key const& key);
PlanF* get_vendor_plan_f (Key const& key);

void add_vendor_plan_d (Key const& key, PlanD plan);
void add_vendor_plan_f (Key const& key, PlanF plan);

template <typename T>
template <Direction D, int M>
void Plan<T>::init_r2c (IntVectND<M> const& fft_size, void* pbf, void* pbb, bool cache)
{
    static_assert(D == Direction::forward || D == Direction::backward);

    kind = (D == Direction::forward) ? Kind::r2c_f : Kind::r2c_b;
    defined = true;
    pf = pbf;
    pb = pbb;

    n = 1;
    for (auto s : fft_size) { n *= s; }
    howmany = 1;

#if defined(AMREX_USE_GPU)
    Key key = {fft_size.template expand<3>(), D, kind};
    if (cache) {
        VendorPlan* cached_plan = nullptr;
        if constexpr (std::is_same_v<float,T>) {
            cached_plan = get_vendor_plan_f(key);
        } else {
            cached_plan = get_vendor_plan_d(key);
        }
        if (cached_plan) {
            plan = *cached_plan;
            return;
        }
    }
#else
    amrex::ignore_unused(cache);
#endif

#if defined(AMREX_USE_CUDA)

    AMREX_CUFFT_SAFE_CALL(cufftCreate(&plan));
    AMREX_CUFFT_SAFE_CALL(cufftSetAutoAllocation(plan, 0));
    cufftType type;
    if constexpr (D == Direction::forward) {
        type = std::is_same_v<float,T> ? CUFFT_R2C : CUFFT_D2Z;
    } else {
        type = std::is_same_v<float,T> ? CUFFT_C2R : CUFFT_Z2D;
    }
    std::size_t work_size;
    if constexpr (M == 1) {
        AMREX_CUFFT_SAFE_CALL
            (cufftMakePlan1d(plan, fft_size[0], type, howmany, &work_size));
    } else if constexpr (M == 2) {
        AMREX_CUFFT_SAFE_CALL
            (cufftMakePlan2d(plan, fft_size[1], fft_size[0], type, &work_size));
    } else if constexpr (M == 3) {
        AMREX_CUFFT_SAFE_CALL
            (cufftMakePlan3d(plan, fft_size[2], fft_size[1], fft_size[0], type, &work_size));
    }

#elif defined(AMREX_USE_HIP)

    auto prec = std::is_same_v<float,T> ? rocfft_precision_single : rocfft_precision_double;
    std::size_t length[M];
    for (int idim = 0; idim < M; ++idim) { length[idim] = fft_size[idim]; }
    if constexpr (D == Direction::forward) {
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&plan, rocfft_placement_notinplace,
                                rocfft_transform_type_real_forward, prec, M,
                                length, howmany, nullptr));
    } else {
        AMREX_ROCFFT_SAFE_CALL
            (rocfft_plan_create(&plan, rocfft_placement_notinplace,
                                rocfft_transform_type_real_inverse, prec, M,
                                length, howmany, nullptr));
    }

#elif defined(AMREX_USE_SYCL)

    mkl_desc_r* pp;
    if (M == 1) {
        pp = new mkl_desc_r(fft_size[0]);
    } else {
        std::vector<std::int64_t> len(M);
        for (int idim = 0; idim < M; ++idim) {
            len[idim] = fft_size[M-1-idim];
        }
        pp = new mkl_desc_r(len);
    }
#ifndef AMREX_USE_MKL_DFTI_2024
    pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT,
                  oneapi::mkl::dft::config_value::NOT_INPLACE);
#else
    pp->set_value(oneapi::mkl::dft::config_param::PLACEMENT, DFTI_NOT_INPLACE);
#endif

    std::vector<std::int64_t> strides(M+1);
    strides[0] = 0;
    strides[M] = 1;
    for (int i = M-1; i >= 1; --i) {
        strides[i] = strides[i+1] * fft_size[M-1-i];
    }

#ifndef AMREX_USE_MKL_DFTI_2024
    pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides);
    pp->set_value(oneapi::mkl::dft::config_param::BWD_STRIDES, strides);
#else
    pp->set_value(oneapi::mkl::dft::config_param::FWD_STRIDES, strides.data());
    // Do not set BWD_STRIDES
#endif
    pp->set_value(oneapi::mkl::dft::config_param::WORKSPACE,
                  oneapi::mkl::dft::config_value::WORKSPACE_EXTERNAL);
    pp->commit(amrex::Gpu::Device::streamQueue());
    plan = pp;

#else /* FFTW */

    if (pf == nullptr || pb == nullptr) {
        defined = false;
        return;
    }

    int size_for_row_major[M];
    for (int idim = 0; idim < M; ++idim) {
        size_for_row_major[idim] = fft_size[M-1-idim];
    }

    if constexpr (std::is_same_v<float,T>) {
        if constexpr (D == Direction::forward) {
            plan = fftwf_plan_dft_r2c
                (M, size_for_row_major, (float*)pf, (fftwf_complex*)pb,
                 FFTW_ESTIMATE);
        } else {
            plan = fftwf_plan_dft_c2r
                (M, size_for_row_major, (fftwf_complex*)pb, (float*)pf,
                 FFTW_ESTIMATE);
        }
    } else {
        if constexpr (D == Direction::forward) {
            plan = fftw_plan_dft_r2c
                (M, size_for_row_major, (double*)pf, (fftw_complex*)pb,
                 FFTW_ESTIMATE);
        } else {
            plan = fftw_plan_dft_c2r
                (M, size_for_row_major, (fftw_complex*)pb, (double*)pf,
                 FFTW_ESTIMATE);
        }
    }
#endif

#if defined(AMREX_USE_GPU)
    if (cache) {
        if constexpr (std::is_same_v<float,T>) {
            add_vendor_plan_f(key, plan);
        } else {
            add_vendor_plan_d(key, plan);
        }
    }
#endif
}

namespace detail
{
    DistributionMapping make_iota_distromap (Long n);

    template <typename FA>
    typename FA::FABType::value_type * get_fab (FA& fa)
    {
        auto myproc = ParallelContext::MyProcSub();
        if (myproc < fa.size()) {
            return fa.fabPtr(myproc);
        } else {
            return nullptr;
        }
    }

    template <typename FA1, typename FA2>
    std::unique_ptr<char,DataDeleter> make_mfs_share (FA1& fa1, FA2& fa2)
    {
        bool not_same_fa = true;
        if constexpr (std::is_same_v<FA1,FA2>) {
            not_same_fa = (&fa1 != &fa2);
        }
        using FAB1 = typename FA1::FABType::value_type;
        using FAB2 = typename FA2::FABType::value_type;
        using T1 = typename FAB1::value_type;
        using T2 = typename FAB2::value_type;
        auto myproc = ParallelContext::MyProcSub();
        bool alloc_1 = (myproc < fa1.size());
        bool alloc_2 = (myproc < fa2.size()) && not_same_fa;
        void* p = nullptr;
        if (alloc_1 && alloc_2) {
            Box const& box1 = fa1.fabbox(myproc);
            Box const& box2 = fa2.fabbox(myproc);
            int ncomp1 = fa1.nComp();
            int ncomp2 = fa2.nComp();
            p = The_Arena()->alloc(std::max(sizeof(T1)*box1.numPts()*ncomp1,
                                            sizeof(T2)*box2.numPts()*ncomp2));
            fa1.setFab(myproc, FAB1(box1, ncomp1, (T1*)p));
            fa2.setFab(myproc, FAB2(box2, ncomp2, (T2*)p));
        } else if (alloc_1) {
            Box const& box1 = fa1.fabbox(myproc);
            int ncomp1 = fa1.nComp();
            p = The_Arena()->alloc(sizeof(T1)*box1.numPts()*ncomp1);
            fa1.setFab(myproc, FAB1(box1, ncomp1, (T1*)p));
        } else if (alloc_2) {
            Box const& box2 = fa2.fabbox(myproc);
            int ncomp2 = fa2.nComp();
            p = The_Arena()->alloc(sizeof(T2)*box2.numPts()*ncomp2);
            fa2.setFab(myproc, FAB2(box2, ncomp2, (T2*)p));
        } else {
            return nullptr;
        }
        return std::unique_ptr<char,DataDeleter>((char*)p, DataDeleter{The_Arena()});
    }
}

struct Swap01
{
    [[nodiscard]] constexpr Dim3 operator() (Dim3 i) const noexcept
    {
        return {i.y, i.x, i.z};
    }

    static constexpr Dim3 Inverse (Dim3 i)
    {
        return {i.y, i.x, i.z};
    }

    [[nodiscard]] constexpr IndexType operator() (IndexType it) const noexcept
    {
        return it;
    }

    static constexpr IndexType Inverse (IndexType it)
    {
        return it;
    }
};

struct Swap02
{
    [[nodiscard]] constexpr Dim3 operator() (Dim3 i) const noexcept
    {
        return {i.z, i.y, i.x};
    }

    static constexpr Dim3 Inverse (Dim3 i)
    {
        return {i.z, i.y, i.x};
    }

    [[nodiscard]] constexpr IndexType operator() (IndexType it) const noexcept
    {
        return it;
    }

    static constexpr IndexType Inverse (IndexType it)
    {
        return it;
    }
};

struct RotateFwd
{
    // dest -> src: (x,y,z) -> (y,z,x)
    [[nodiscard]] constexpr Dim3 operator() (Dim3 i) const noexcept
    {
        return {i.y, i.z, i.x};
    }

    // src -> dest: (x,y,z) -> (z,x,y)
    static constexpr Dim3 Inverse (Dim3 i)
    {
        return {i.z, i.x, i.y};
    }

    [[nodiscard]] constexpr IndexType operator() (IndexType it) const noexcept
    {
        return it;
    }

    static constexpr IndexType Inverse (IndexType it)
    {
        return it;
    }
};

struct RotateBwd
{
    // dest -> src: (x,y,z) -> (z,x,y)
    [[nodiscard]] constexpr Dim3 operator() (Dim3 i) const noexcept
    {
        return {i.z, i.x, i.y};
    }

    // src -> dest: (x,y,z) -> (y,z,x)
    static constexpr Dim3 Inverse (Dim3 i)
    {
        return {i.y, i.z, i.x};
    }

    [[nodiscard]] constexpr IndexType operator() (IndexType it) const noexcept
    {
        return it;
    }

    static constexpr IndexType Inverse (IndexType it)
    {
        return it;
    }
};

}

#endif
