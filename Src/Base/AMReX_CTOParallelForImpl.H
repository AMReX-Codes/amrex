#ifndef AMREX_CTO_PARALLEL_FOR_H_
#define AMREX_CTO_PARALLEL_FOR_H_

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_TypeList.H>

#include <array>
#include <type_traits>

/* This header is not for the users to include directly.  It's meant to be
 * included in AMReX_GpuLaunch.H, which has included the headers needed
 * here. */

/* Thank Maikel Nadolski and Alex Sinn for the techniques used here! */

namespace amrex {

template <int... ctr>
struct CompileTimeOptions {
    // TypeList is defined in AMReX_TypeList.H
    using list_type = TypeList<std::integral_constant<int, ctr>...>;
};

#if (__cplusplus >= 201703L)

namespace detail
{
    template<class F, int... ctr>
    struct CTOWrapper {
        F f;

        template<class... Args>
        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        auto operator() (Args... args) const noexcept
            -> decltype(f(args..., std::integral_constant<int, ctr>{}...)) {
            return f(args..., std::integral_constant<int, ctr>{}...);
        }

        AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
        static constexpr
        std::array<int, sizeof...(ctr)> GetOptions () noexcept {
            return {ctr...};
        }
    };

    template <class L, class F, typename... As>
    bool
    AnyCTO_helper2 (const L& l, const F& f, TypeList<As...>,
                    std::array<int,sizeof...(As)> const& runtime_options)
    {
        if (runtime_options == std::array<int,sizeof...(As)>{As::value...}) {
            l(CTOWrapper<F, As::value...>{f});
            return true;
        } else {
            return false;
        }
    }

    template <class L, class F, typename... PPs, typename RO>
    void
    AnyCTO_helper1 (const L& l, const F& f, TypeList<PPs...>, RO const& runtime_options)
    {
        bool found_option = (false || ... || AnyCTO_helper2(l, f, PPs{}, runtime_options));
        amrex::ignore_unused(found_option);
        AMREX_ASSERT(found_option);
    }
}

#endif

/**
 * \brief Compile time optimization of kernels with run time options.
 *
 * This is a generalized version of ParallelFor with CTOs that can support any function that
 * takes in one lambda to launch a GPU kernel such as ParallelFor, ParallelForRNG, launch, etc.
 * It uses fold expression to generate kernel launches for all combinations
 * of the run time options.  The kernel function can use constexpr if to
 * discard unused code blocks for better run time performance.  In the
 * example below, the code will be expanded into 4*2=8 normal ParallelForRNGs
 * for all combinations of the run time parameters.
 \verbatim
    int A_runtime_option = ...;
    int B_runtime_option = ...;
    enum A_options : int { A0, A1, A2, A3 };
    enum B_options : int { B0, B1 };
    AnyCTO(TypeList<CompileTimeOptions<A0,A1,A2,A3>,
                    CompileTimeOptions<B0,B1>>{},
        {A_runtime_option, B_runtime_option},
        [&](auto cto_func){
            ParallelForRNG(N, cto_func);
        },
        [=] AMREX_GPU_DEVICE (int i, const RandomEngine& engine,
                              auto A_control, auto B_control)
        {
            ...
            if constexpr (A_control.value == A0) {
                ...
            } else if constexpr (A_control.value == A1) {
                ...
            } else if constexpr (A_control.value == A2) {
                ...
            } else {
                ...
            }
            if constexpr (A_control.value != A3 && B_control.value == B1) {
                ...
            }
            ...
        }
    );

    constexpr int nthreads_per_block = ...;
    int nblocks = ...;
    AnyCTO(TypeList<CompileTimeOptions<A0,A1,A2,A3>,
                    CompileTimeOptions<B0,B1>>{},
        {A_runtime_option, B_runtime_option},
        [&](auto cto_func){
            launch<nthreads_per_block>(nblocks, Gpu::gpuStream(), cto_func);
        },
        [=] AMREX_GPU_DEVICE (auto A_control, auto B_control){
            ...
        }
    );
 \endverbatim
 * The static member function cto_func.GetOptions() can be used to obtain the runtime_options
 * passed into AnyCTO, but at compile time. This enables some advanced use cases,
 * such as changing the number of threads per block or the dimensionality of ParallelFor at runtime.
 * For the second example -> decltype(void(intvect.size())) is necessary to
 * disambiguate IntVectND<1> and int for the first argument of the kernel function.
 \verbatim
    int nthreads_per_block = ...;
    AnyCTO(TypeList<CompileTimeOptions<128,256,512,1024>>{},
        {nthreads_per_block},
        [&](auto cto_func){
            constexpr std::array<int, 1> ctos = cto_func.GetOptions();
            constexpr int c_nthreads_per_block = ctos[0];
            ParallelFor<c_nthreads_per_block>(N, cto_func);
        },
        [=] AMREX_GPU_DEVICE (int i, auto){
            ...
        }
    );

    BoxND<6> box6D = ...;
    int dims_needed = ...;
    AnyCTO(TypeList<CompileTimeOptions<1,2,3,4,5,6>>{},
        {dims_needed},
        [&](auto cto_func){
            constexpr std::array<int, 1> ctos = cto_func.GetOptions();
            constexpr int c_dims_needed = ctos[0];
            const auto box = BoxShrink<c_dims_needed>(box6D);
            ParallelFor(box, cto_func);
        },
        [=] AMREX_GPU_DEVICE (auto intvect, auto) -> decltype(void(intvect.size())) {
            ...
        }
    );
 \endverbatim

 * Note that due to a limitation of CUDA's extended device lambda, the
 * constexpr if block cannot be the one that captures a variable first.
 * If nvcc complains about it, you will have to manually capture it outside
 * constexpr if. Alternatively, the constexpr if can be replaced with a regular if.
 * Compilers can still perform the same optimizations since the condition is known at compile time.
 * The data type for the parameters is int.
 *
 * \param list_of_compile_time_options list of all possible values of the parameters.
 * \param runtime_options the run time parameters.
 * \param l a callable object containing a CPU function that launches the provided GPU kernel.
 * \param f a callable object containing the GPU kernel with optimizations.
 */
template <class L, class F, typename... CTOs>
void AnyCTO ([[maybe_unused]] TypeList<CTOs...> list_of_compile_time_options,
             std::array<int,sizeof...(CTOs)> const& runtime_options,
             L&& l, F&& f)
{
#if (__cplusplus >= 201703L)
    detail::AnyCTO_helper1(std::forward<L>(l), std::forward<F>(f),
                           CartesianProduct(typename CTOs::list_type{}...),
                           runtime_options);
#else
    amrex::ignore_unused(runtime_options, l, f);
    static_assert(std::is_integral<F>::value, "This requires C++17");
#endif
}

template <int MT, typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral_v<T>>
ParallelFor (TypeList<CTOs...> ctos,
             std::array<int,sizeof...(CTOs)> const& runtime_options,
             T N, F&& f)
{
    AnyCTO(ctos, runtime_options,
        [&](auto cto_func){
            ParallelFor<MT>(N, cto_func);
        },
        std::forward<F>(f)
    );
}

template <int MT, class F, int dim, typename... CTOs>
void ParallelFor (TypeList<CTOs...> ctos,
                  std::array<int,sizeof...(CTOs)> const& runtime_options,
                  BoxND<dim> const& box, F&& f)
{
    AnyCTO(ctos, runtime_options,
        [&](auto cto_func){
            ParallelFor<MT>(box, cto_func);
        },
        std::forward<F>(f)
    );
}

template <int MT, typename T, class F, int dim, typename... CTOs>
std::enable_if_t<std::is_integral_v<T>>
ParallelFor (TypeList<CTOs...> ctos,
             std::array<int,sizeof...(CTOs)> const& runtime_options,
             BoxND<dim> const& box, T ncomp, F&& f)
{
    AnyCTO(ctos, runtime_options,
        [&](auto cto_func){
            ParallelFor<MT>(box, ncomp, cto_func);
        },
        std::forward<F>(f)
    );
}

/**
 * \brief ParallelFor with compile time optimization of kernels with run time options.
 *
 * It uses fold expression to generate kernel launches for all combinations
 * of the run time options.  The kernel function can use constexpr if to
 * discard unused code blocks for better run time performance.  In the
 * example below, the code will be expanded into 4*2=8 normal ParallelFors
 * for all combinations of the run time parameters.
 \verbatim
     int A_runtime_option = ...;
     int B_runtime_option = ...;
     enum A_options : int { A0, A1, A2, A3};
     enum B_options : int { B0, B1 };
     ParallelFor(TypeList<CompileTimeOptions<A0,A1,A2,A3>,
                          CompileTimeOptions<B0,B1>>{},
                 {A_runtime_option, B_runtime_option},
                 N, [=] AMREX_GPU_DEVICE (int i, auto A_control, auto B_control)
     {
         ...
         if constexpr (A_control.value == A0) {
             ...
         } else if constexpr (A_control.value == A1) {
             ...
         } else if constexpr (A_control.value == A2) {
             ...
         } else {
             ...
         }
         if constexpr (A_control.value != A3 && B_control.value == B1) {
             ...
         }
         ...
     });
 \endverbatim
 * Note that due to a limitation of CUDA's extended device lambda, the
 * constexpr if block cannot be the one that captures a variable first.
 * If nvcc complains about it, you will have to manually capture it outside
 * constexpr if.  The data type for the parameters is int.
 *
 * \param ctos   list of all possible values of the parameters.
 * \param option the run time parameters.
 * \param N      an integer specifying the 1D for loop's range.
 * \param f      a callable object taking an integer and working on that iteration.
 */
template <typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral_v<T>>
ParallelFor (TypeList<CTOs...> ctos,
             std::array<int,sizeof...(CTOs)> const& option,
             T N, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(ctos, option, N, std::forward<F>(f));
}

/**
 * \brief ParallelFor with compile time optimization of kernels with run time options.
 *
 * It uses fold expression to generate kernel launches for all combinations
 * of the run time options.  The kernel function can use constexpr if to
 * discard unused code blocks for better run time performance.  In the
 * example below, the code will be expanded into 4*2=8 normal ParallelFors
 * for all combinations of the run time parameters.
 \verbatim
     int A_runtime_option = ...;
     int B_runtime_option = ...;
     enum A_options : int { A0, A1, A2, A3};
     enum B_options : int { B0, B1 };
     ParallelFor(TypeList<CompileTimeOptions<A0,A1,A2,A3>,
                          CompileTimeOptions<B0,B1>>{},
                 {A_runtime_option, B_runtime_option},
                 box, [=] AMREX_GPU_DEVICE (int i, int j, int k,
                                            auto A_control, auto B_control)
     {
         ...
         if constexpr (A_control.value == A0) {
             ...
         } else if constexpr (A_control.value == A1) {
             ...
         } else if constexpr (A_control.value == A2) {
             ...
         } else {
             ...
         }
         if constexpr (A_control.value != A3 && B_control.value == B1) {
             ...
         }
         ...
     });
 \endverbatim
 * Note that due to a limitation of CUDA's extended device lambda, the
 * constexpr if block cannot be the one that captures a variable first.
 * If nvcc complains about it, you will have to manually capture it outside
 * constexpr if.  The data type for the parameters is int.
 *
 * \param ctos   list of all possible values of the parameters.
 * \param option the run time parameters.
 * \param box      a Box specifying the 3D for loop's range.
 * \param f        a callable object taking three integers and working on the given cell.
 */
template <class F, int dim, typename... CTOs>
void ParallelFor (TypeList<CTOs...> ctos,
                  std::array<int,sizeof...(CTOs)> const& option,
                  BoxND<dim> const& box, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(ctos, option, box, std::forward<F>(f));
}

/**
 * \brief ParallelFor with compile time optimization of kernels with run time options.
 *
 * It uses fold expression to generate kernel launches for all combinations
 * of the run time options.  The kernel function can use constexpr if to
 * discard unused code blocks for better run time performance.  In the
 * example below, the code will be expanded into 4*2=8 normal ParallelFors
 * for all combinations of the run time parameters.
 \verbatim
     int A_runtime_option = ...;
     int B_runtime_option = ...;
     enum A_options : int { A0, A1, A2, A3};
     enum B_options : int { B0, B1 };
     ParallelFor(TypeList<CompileTimeOptions<A0,A1,A2,A3>,
                          CompileTimeOptions<B0,B1>>{},
                 {A_runtime_option, B_runtime_option},
                 box, ncomp, [=] AMREX_GPU_DEVICE (int i, int j, int k, int n,
                                                   auto A_control, auto B_control)
     {
         ...
         if constexpr (A_control.value == A0) {
             ...
         } else if constexpr (A_control.value == A1) {
             ...
         } else if constexpr (A_control.value == A2) {
             ...
         } else {
             ...
         }
         if constexpr (A_control.value != A3 && B_control.value == B1) {
             ...
         }
         ...
     });
 \endverbatim
 * Note that due to a limitation of CUDA's extended device lambda, the
 * constexpr if block cannot be the one that captures a variable first.
 * If nvcc complains about it, you will have to manually capture it outside
 * constexpr if.  The data type for the parameters is int.
 *
 * \param ctos   list of all possible values of the parameters.
 * \param option the run time parameters.
 * \param box    a Box specifying the iteration in 3D space.
 * \param ncomp  an integer specifying the range for iteration over components.
 * \param f      a callable object taking three integers and working on the given cell.
 */
template <typename T, class F, int dim, typename... CTOs>
std::enable_if_t<std::is_integral_v<T>>
ParallelFor (TypeList<CTOs...> ctos,
             std::array<int,sizeof...(CTOs)> const& option,
             BoxND<dim> const& box, T ncomp, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(ctos, option, box, ncomp, std::forward<F>(f));
}

}

#endif
