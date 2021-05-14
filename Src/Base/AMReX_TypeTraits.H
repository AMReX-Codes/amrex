#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_
#include <AMReX_Config.H>

#include <AMReX_Extension.H>
#include <type_traits>

// In case they are still used by applications
#define AMREX_IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable<T>::value
#define AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T) std::is_trivially_default_constructible<T>::value

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
    struct MaybeHostDeviceRunnable<T, std::enable_if_t<__nv_is_extended_device_lambda_closure_type(T)> >
        : std::false_type {};

    template <class T>
    struct DefinitelyNotHostRunnable<T, std::enable_if_t<__nv_is_extended_device_lambda_closure_type(T)> >
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
    //                                                           [traits.IsDetected]
    //
    // We use IsDetected as a SFINAE tool to test for valid expressions.
    //
    // is_detected was proposed to C++ but is surpassed by concepts.
    //
    // The implementation is taken from
    //
    // https://en.cppreference.com/w/cpp/experimental/is_detected

    namespace detail {
        template <class...> using Void_t = void;

        struct Nonesuch {
            Nonesuch() = delete;
            ~Nonesuch() = delete;
            Nonesuch(Nonesuch const&) = delete;
            void operator=(Nonesuch const&) = delete;
        };

        template <class Default, class AlwaysVoid, template <class...> class Op,
                class... Args>
        struct Detector {
            using value_t = std::false_type;
            using type = Default;
        };

        template <class Default, template <class...> class Op, class... Args>
        struct Detector<Default, Void_t<Op<Args...>>, Op, Args...> {
            using value_t = std::true_type;
            using type = Op<Args...>;
        };
    }

    template <template <class...> class Op, class... Args>
    using IsDetected = typename detail::Detector<detail::Nonesuch, void, Op, Args...>::value_t;

    template <template <class...> class Op, class... Args>
    using Detected_t = typename detail::Detector<detail::Nonesuch, void, Op, Args...>::type;

    template <class Default, template <class...> class Op, class... Args>
    using DetectedOr = typename detail::Detector<Default, void, Op, Args...>::type;

    template <class Expected, template <typename...> class Op, class... Args>
    using IsDetectedExact = std::is_same<Expected, Detected_t<Op, Args...>>;

    ////////////////////////////////////////////////////////////////////////////////
    //                                                           [traits.IsCallable]

    namespace detail {
        template <typename T, typename... Args>
        using call_result_t = decltype(std::declval<T>()(std::declval<Args>()...));
    }

    //! \brief Test if a given type T is callable with arguments of type Args...
    //!
    //! This type trait is different from std::is_invocable since it only cares for the call syntax
    //! f(args...) and disregards pointers to class methods and such.
    //!
    //! For consistency we use IsCallable over is_invocable since AMReX relies on call syntax
    //! everywhere else. Using an equivalent function to std::invoke would increase call stacks in
    //! debug mode, add more template instantiations and will obfuscate lots of generic code that
    //! takes functions.
    template <typename T, typename... Args>
    struct IsCallable : IsDetected<detail::call_result_t, T, Args...> {};

    //! \brief Test if a given type T is callable with arguments of type Args...
    //!
    //! \see IsCallable
    template <typename R, typename T, typename... Args>
    struct IsCallableR : IsDetectedExact<R, detail::call_result_t, T, Args...> {};

    ////////////////////////////////////////////////////////////////////////////////
    //                                                          [traits.Conjunction]
    //                                                          [traits.Disjunction]
    //                                                          [traits.Negation]


    //! Logical traits let us combine multiple type requirements in one enable_if_t clause.
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

    template <class B>
    using Negation = std::integral_constant<bool, !bool(B::value)>;
#endif

    // Move this down, because doxygen can not parse anything below IsStoreAtomic
    template <class T, class Enable = void>
    struct IsStoreAtomic : std::false_type {};
    //
    template <class T>
    struct IsStoreAtomic<T, typename std::enable_if <
                                std::is_arithmetic<T>::value
                                && sizeof(T) <= 8 >::type >
        : std::true_type {};

}

#endif
