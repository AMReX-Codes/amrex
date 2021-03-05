#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_
#include <AMReX_Config.H>

#include <type_traits>
#include <AMReX_Extension.H>

// workaround missing "is_trivially_copyable" in g++ < 5.0
#if defined(__GNUG__) && __GNUC__ < 5 && !defined(AMREX_CXX_CLANG)
#define AMREX_IS_TRIVIALLY_COPYABLE(T) __has_trivial_copy(T)
#define AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T) std::has_trivial_default_constructor<T>::value
#else
#define AMREX_IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#define AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T) std::is_trivially_default_constructible<T>::value
#endif

namespace amrex
{
    template <class T> class BaseFab;
    template <class FAB> class FabArray;

    template <class A, class Enable = void> struct IsBaseFab : std::false_type {};
    //
    template <class D>
    struct IsBaseFab<D, typename std::enable_if<
                            std::is_base_of<BaseFab<typename D::value_type>,
                                            D>::value>::type>
        : std::true_type {};


    template <class A, class Enable = void> struct IsFabArray : std::false_type {};
    //
    template <class D>
    struct IsFabArray<D, typename std::enable_if<
                             std::is_base_of<FabArray<typename D::FABType::value_type>,
                                             D>::value>::type>
        : std::true_type {};

    template <bool B, class T = void>
    using EnableIf_t = typename std::enable_if<B,T>::type;


    template <class T, class Enable = void>
    struct IsStoreAtomic : std::false_type {};
    //
    template <class T>
    struct IsStoreAtomic<T, typename std::enable_if <
                                std::is_arithmetic<T>::value
                                && sizeof(T) <= 8 >::type >
        : std::true_type {};

    template <class T, class Enable = void>
    struct HasAtomicAdd : std::false_type {};
    template <> struct HasAtomicAdd<int> : std::true_type {};
    template <> struct HasAtomicAdd<long> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned int> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned long long> : std::true_type {};
    template <> struct HasAtomicAdd<float> : std::true_type {};
    template <> struct HasAtomicAdd<double> : std::true_type {};

#ifdef AMREX_PARTICLES
    template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
              template<class> class Allocator>
    class ParIterBase;

    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
              template<class> class Allocator>
    class ParIter;

    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
              template<class> class Allocator>
    class ParConstIter;

    class MFIter;
    class ParticleContainerBase;

    template <typename T>
    struct IsParticleIterator : public std::is_base_of<MFIter, T>::type {};

    template <typename T>
    struct IsParticleContainer : public std::is_base_of<ParticleContainerBase, T>::type {};
#endif

#ifdef AMREX_USE_GPU

    template <class T, class Enable = void>
    struct MaybeDeviceRunnable : std::true_type {};

    template <class T, class Enable = void>
    struct MaybeHostDeviceRunnable : std::true_type {};

    template <class T, class Enable = void>
    struct DefinitelyNotHostRunnable : std::false_type {};

#if defined(AMREX_USE_CUDA)

    template <class T>
    struct MaybeHostDeviceRunnable<T, amrex::EnableIf_t<__nv_is_extended_device_lambda_closure_type(T)> >
        : std::false_type {};

    template <class T>
    struct DefinitelyNotHostRunnable<T, amrex::EnableIf_t<__nv_is_extended_device_lambda_closure_type(T)> >
        : std::true_type {};

#elif defined(AMREX_USE_HIP)

    // xxxxx HIP todo

#endif

#endif

    template <typename T, typename U1, typename... Us>
    struct Same;

    template <typename T, typename U>
    struct Same<T,U>
    {
        static constexpr bool value = std::is_same<T,U>::value;
    };

    template <typename T, typename U1, typename... Us>
    struct Same
    {
        static constexpr bool value = std::is_same<T,U1>::value && Same<T,Us...>::value;
    };

    ////////////////////////////////////////////////////////////////////////////////
    //                                                          [traits.is_detected]

    /// \defgroup type-traits Type Traits
    /// Type traits are used to enforce requirements on types with SFINAE.

    template <class...> using void_t = void;

    template <typename T> struct nodeduce { using type = T; };
    template <typename T> using nodeduce_t = typename nodeduce<T>::type;

    struct nonesuch {
        nonesuch() = delete;
        ~nonesuch() = delete;
        nonesuch(nonesuch const&) = delete;
        void operator=(nonesuch const&) = delete;
    };

    namespace detail {
        template <class Default, class AlwaysVoid, template <class...> class Op,
                class... Args>
        struct detector {
            using value_t = std::false_type;
            using type = Default;
        };

        template <class Default, template <class...> class Op, class... Args>
        struct detector<Default, void_t<Op<Args...>>, Op, Args...> {
            using value_t = std::true_type;
            using type = Op<Args...>;
        };

        template <template <class...> class Op, class... Args>
        using is_detected =
            typename detail::detector<nonesuch, void, Op, Args...>::value_t;

        /// \ingroup type-traits
        /// Returns the type of `Op<Args...>` or `nonesuch`
        template <template <class...> class Op, class... Args>
        using detected_t = typename detail::detector<nonesuch, void, Op, Args...>::type;

        template <class Default, template <class...> class Op, class... Args>
        using detected_or = detail::detector<Default, void, Op, Args...>;

        template <class Default, template <class...> class Op, class... Args>
        using detected_or_t = typename detected_or<Default, Op, Args...>::type;

        template <class Expected, template <typename...> class Op, class... Args>
        using is_detected_exact = std::is_same<Expected, detected_t<Op, Args...>>;
    } // namespace detail

    /// \ingroup type-traits
    /// This is `std::true_type` if `Op<Args...>` is a valid SFINAE expression.
    template <template <class...> class Op, class... Args>
    struct IsDetected : detail::is_detected<Op, Args...> {};

    /// \ingroup type-traits
    /// Returns the type of `Op<Args...>` or `nonesuch`
    template <template <class...> class Op, class... Args>
    using Detected_t = detail::detected_t<Op, Args...>;

    /// \ingroup type-traits
    /// Returns the type of `Op<Args...>` or `Default`
    template <class Default, template <class...> class Op, class... Args>
    struct DetectedOr : detail::detected_or<Default, Op, Args...> {};

    /// \ingroup type-traits
    /// This is `std::true_type` if `Op<Args...>` is a valid SFINAE expression and
    /// the return type is exactly `Expected`.
    template <class Expected, template <typename...> class Op, class... Args>
    struct IsDetectedExact : detail::is_detected_exact<Expected, Op, Args...> {};

    namespace detail {
        template <typename T, typename... Args>
        using call_result_t = decltype(std::declval<T>()(std::declval<Args>()...));
    }

    /// @{
    /// \ingroup type-traits
    /// This is different from std::is_invocable since it only cares for the call syntax f(args...)
    /// and disregards pointers to class methods and such.
    template <typename T, typename... Args>
    struct IsCallable : IsDetected<detail::call_result_t, T, Args...> {};

    template <typename R, typename T, typename... Args>
    struct IsCallableR : IsDetectedExact<R, detail::call_result_t, T, Args...> {};
    /// @}

    #if defined(__cpp_lib_logical_traits)
    template <typename... Args> using Conjunction = std::conjunction<Args...>;
    template <typename... Args> using Disjunction = std::disjunction<Args...>;
    template <typename... Args> using Negation = std::negation<Args...>;
    #elif defined(__cpp_lib_experimental_logical_traits)
    template <typename... Args> using Conjunction = std::experimental::conjunction<Args...>;
    template <typename... Args> using Disjunction = std::experimental::disjunction<Args...>;
    template <typename... Args> using Negation = std::experimental::negation<Args...>;
    #else
    template <class...> struct Conjunction : std::true_type {};
    template <class B1> struct Conjunction<B1> : B1 {};
    template <class B1, class... Bn>
    struct Conjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), Conjunction<Bn...>, B1> {};

    template <class...> struct Disjunction : std::false_type {};
    template <class B1> struct Disjunction<B1> : B1 {};
    template <class B1, class... Bn>
    struct Disjunction<B1, Bn...>
    : std::conditional_t<bool(B1::value), B1, Disjunction<Bn...>> {};

    template <class Bool>
    using Negation = std::integral_constant<bool, !bool(Bool::value)>;
    #endif
}

#endif
