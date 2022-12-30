#ifndef AMREX_CTO_PARALLEL_FOR_H_
#define AMREX_CTO_PARALLEL_FOR_H_

#include <AMReX_BLassert.H>
#include <AMReX_Box.H>
#include <AMReX_Tuple.H>

#include <array>
#include <type_traits>

/* This header is not for the users to include directly.  It's meant to be
 * included in AMReX_GpuLaunch.H, which has included the headers needed
 * here. */

/* Thank Maikel Nadolski and Alex Sinn for the techniques used here! */

namespace amrex {

template <int... ctr>
struct CompileTimeOptions {
    // TypeList is defined in AMReX_Tuple.H
    using list_type = TypeList<std::integral_constant<int, ctr>...>;
};

#if (__cplusplus >= 201703L)

namespace detail
{
    template <int MT, typename T, class F, typename... As>
    std::enable_if_t<std::is_integral<T>::value || std::is_same<T,Box>::value, bool>
    ParallelFor_helper2 (T const& N, F&& f, TypeList<As...>,
                         std::array<int,sizeof...(As)> const& runtime_options)
    {
        if (runtime_options == std::array<int,sizeof...(As)>{As::value...}) {
            if constexpr (std::is_integral<T>::value) {
                ParallelFor<MT>(N, [f] AMREX_GPU_DEVICE (T i) noexcept
                {
                    f(i, As{}...);
                });
            } else {
                ParallelFor<MT>(N, [f] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
                {
                    f(i, j, k, As{}...);
                });
            }
            return true;
        } else {
            return false;
        }
    }

    template <int MT, typename T, class F, typename... As>
    std::enable_if_t<std::is_integral<T>::value, bool>
    ParallelFor_helper2 (Box const& box, T ncomp, F&& f, TypeList<As...>,
                         std::array<int,sizeof...(As)> const& runtime_options)
    {
        if (runtime_options == std::array<int,sizeof...(As)>{As::value...}) {
            ParallelFor<MT>(box, ncomp, [f] AMREX_GPU_DEVICE (int i, int j, int k, T n) noexcept
            {
                f(i, j, k, n, As{}...);
            });
            return true;
        } else {
            return false;
        }
    }

    template <int MT, typename T, class F, typename... PPs, typename RO>
    std::enable_if_t<std::is_integral<T>::value || std::is_same<T,Box>::value>
    ParallelFor_helper1 (T const& N, F&& f, TypeList<PPs...>,
                         RO const& runtime_options)
    {
        bool found_option = (false || ... ||
                             ParallelFor_helper2<MT>(N, std::forward<F>(f),
                                                     PPs{}, runtime_options));
        amrex::ignore_unused(found_option);
        AMREX_ASSERT(found_option);
    }

    template <int MT, typename T, class F, typename... PPs, typename RO>
    std::enable_if_t<std::is_integral<T>::value>
    ParallelFor_helper1 (Box const& box, T ncomp, F&& f, TypeList<PPs...>,
                         RO const& runtime_options)
    {
        bool found_option = (false || ... ||
                             ParallelFor_helper2<MT>(box, ncomp, std::forward<F>(f),
                                                     PPs{}, runtime_options));
        amrex::ignore_unused(found_option);
        AMREX_ASSERT(found_option);
    }
}

#endif

template <int MT, typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral<T>::value>
ParallelFor (TypeList<CTOs...> /*list_of_compile_time_options*/,
             std::array<int,sizeof...(CTOs)> const& runtime_options,
             T N, F&& f)
{
#if (__cplusplus >= 201703L)
    detail::ParallelFor_helper1<MT>(N, std::forward<F>(f),
                                    CartesianProduct(typename CTOs::list_type{}...),
                                    runtime_options);
#else
    amrex::ignore_unused(N, f, runtime_options);
    static_assert(std::is_integral<F>::value, "This requires C++17");
#endif
}

template <int MT, class F, typename... CTOs>
void ParallelFor (TypeList<CTOs...> /*list_of_compile_time_options*/,
                  std::array<int,sizeof...(CTOs)> const& runtime_options,
                  Box const& box, F&& f)
{
#if (__cplusplus >= 201703L)
    detail::ParallelFor_helper1<MT>(box, std::forward<F>(f),
                                    CartesianProduct(typename CTOs::list_type{}...),
                                    runtime_options);
#else
    amrex::ignore_unused(box, f, runtime_options);
    static_assert(std::is_integral<F>::value, "This requires C++17");
#endif
}

template <int MT, typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral<T>::value>
ParallelFor (TypeList<CTOs...> /*list_of_compile_time_options*/,
             std::array<int,sizeof...(CTOs)> const& runtime_options,
             Box const& box, T ncomp, F&& f)
{
#if (__cplusplus >= 201703L)
    detail::ParallelFor_helper1<MT>(box, ncomp, std::forward<F>(f),
                                    CartesianProduct(typename CTOs::list_type{}...),
                                    runtime_options);
#else
    amrex::ignore_unused(box, ncomp, f, runtime_options);
    static_assert(std::is_integral<F>::value, "This requires C++17");
#endif
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
         else {
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
 * \param N      an interger specifying the 1D for loop's range.
 * \param f      a callable object taking an integer and working on that iteration.
 */
template <typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral<T>::value>
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
         else {
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
template <class F, typename... CTOs>
void ParallelFor (TypeList<CTOs...> ctos,
                  std::array<int,sizeof...(CTOs)> const& option,
                  Box const& box, F&& f)
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
         else {
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
template <typename T, class F, typename... CTOs>
std::enable_if_t<std::is_integral<T>::value>
ParallelFor (TypeList<CTOs...> ctos,
             std::array<int,sizeof...(CTOs)> const& option,
             Box const& box, T ncomp, F&& f)
{
    ParallelFor<AMREX_GPU_MAX_THREADS>(ctos, option, box, ncomp, std::forward<F>(f));
}

}

#endif
