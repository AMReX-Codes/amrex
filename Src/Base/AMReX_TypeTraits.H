#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_TypeTraits.H
 * \brief Type traits and metaprogramming utilities used throughout AMReX.
 */

#include <AMReX_Extension.H>
#include <vector>
#include <type_traits>

// In case they are still used by applications
#define AMREX_IS_TRIVIALLY_COPYABLE(T) std::is_trivially_copyable_v<T>
#define AMREX_IS_TRIVIALLY_DEFAULT_CONSTRUCTIBLE(T) std::is_trivially_default_constructible_v<T>

namespace amrex
{
    template <class T> class BaseFab;
    template <class FAB> class FabArray;

    /**
     * \class amrex::IsBaseFab
     * \ingroup amrex_utilities
     * \brief True if D is derived from BaseFab<D::value_type>.
     */
    template <class A, class Enable = void> struct IsBaseFab : std::false_type {};
    //
    template <class D>
    struct IsBaseFab<D, std::enable_if_t<
                            std::is_base_of_v<BaseFab<typename D::value_type>,
                                            D>>>
        : std::true_type {};
    //
    /// \ingroup amrex_utilities
    //! True if A is derived from BaseFab<A::value_type>.
    template <class A>
    inline constexpr bool IsBaseFab_v = IsBaseFab<A>::value;

    /**
     * \class amrex::IsFabArray
     * \ingroup amrex_utilities
     * \brief True if D is derived from FabArray<D::FABType::value_type>.
     */
    template <class A, class Enable = void> struct IsFabArray : std::false_type {};
    //
    template <class D>
    struct IsFabArray<D, std::enable_if_t<
                             std::is_base_of_v<FabArray<typename D::FABType::value_type>,
                                             D>>>
        : std::true_type {};
    //
    /// \ingroup amrex_utilities
    //! True if A is derived from FabArray<A::FABType::value_type>.
    template <class A>
    inline constexpr bool IsFabArray_v = IsFabArray<A>::value;

    /**
     * \class amrex::IsMultiFabLike
     * \ingroup amrex_utilities
     * \brief True if M is both a FabArray and has a BaseFab-derived fab type.
     */
    template <class M, class Enable = void>
    struct IsMultiFabLike : std::false_type {};
    //
    template <class M>
    struct IsMultiFabLike<M, std::enable_if_t<IsFabArray_v<M> &&
                                              IsBaseFab_v<typename M::fab_type> > >
        : std::true_type {};
    //
    /// \ingroup amrex_utilities
    //! True if M is both a FabArray and has a BaseFab-derived fab type.
    template <class M>
    inline constexpr bool IsMultiFabLike_v = IsMultiFabLike<M>::value;


    /// \ingroup amrex_utilities
    //! Alias for std::enable_if_t.
    template <bool B, class T = void>
    using EnableIf_t = std::enable_if_t<B,T>;

    /**
     * \class amrex::HasAtomicAdd
     * \ingroup amrex_utilities
     * \brief True if T supports GPU atomic add operations.
     */
    template <class T, class Enable = void>
    struct HasAtomicAdd : std::false_type {};
    template <> struct HasAtomicAdd<int> : std::true_type {};
    template <> struct HasAtomicAdd<long> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned int> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned long long> : std::true_type {};
    template <> struct HasAtomicAdd<float> : std::true_type {};
    template <> struct HasAtomicAdd<double> : std::true_type {};

    /**
     * \class amrex::IsAddAssignable
     * \ingroup amrex_utilities
     * \brief True if T supports the += operator.
     */
    // support +=?
    template <class T, class Enable = void>
    struct IsAddAssignable : std::false_type {};
    //
    template <class T>
    struct IsAddAssignable<T, std::void_t<decltype(std::declval<T&>() += std::declval<T>())>>
        : std::true_type {};

    class MFIter;
    /**
     * \class amrex::IsMultiFabIterator
     * \ingroup amrex_utilities
     * \brief True if T is derived from MFIter.
     */
    template <typename T>
    struct IsMultiFabIterator : public std::is_base_of<MFIter, T>::type {};

#ifdef AMREX_PARTICLES
    class ParticleContainerBase;

    /**
     * \class amrex::IsParticleIterator
     * \ingroup amrex_utilities
     * \brief True if T is derived from MFIter (used as a proxy for particle iterators).
     */
    template <typename T>
    struct IsParticleIterator : public std::is_base_of<MFIter, T>::type {}; // not exactly right

    /**
     * \class amrex::IsParticleContainer
     * \ingroup amrex_utilities
     * \brief True if T is derived from ParticleContainerBase.
     */
    template <typename T>
    struct IsParticleContainer : public std::is_base_of<ParticleContainerBase, T>::type {};
#endif

    /**
     * \class amrex::DefinitelyNotHostRunnable
     * \ingroup amrex_utilities
     * \brief True if T is definitely not callable on the host
     *        (e.g., a CUDA \c __device__-only lambda).
     */
    template <class T, class Enable = void>
    struct DefinitelyNotHostRunnable : std::false_type {};

#ifdef AMREX_USE_GPU

    /**
     * \class amrex::MaybeDeviceRunnable
     * \ingroup amrex_utilities
     * \brief True if T might be callable inside a GPU device kernel.
     */
    template <class T, class Enable = void>
    struct MaybeDeviceRunnable : std::true_type {};

    /**
     * \class amrex::MaybeHostDeviceRunnable
     * \ingroup amrex_utilities
     * \brief True if T might be callable on both host and device.
     */
    template <class T, class Enable = void>
    struct MaybeHostDeviceRunnable : std::true_type {};

#if defined(AMREX_USE_CUDA) && defined(__NVCC__)

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

    /**
     * \class amrex::Same
     * \ingroup amrex_utilities
     * \brief True if all types U1, Us... are the same as T.
     */
    template <typename T, typename U1, typename... Us>
    struct Same;

    template <typename T, typename U>
    struct Same<T,U>
    {
        static constexpr bool value = std::is_same_v<T,U>;
    };

    template <typename T, typename U1, typename... Us>
    struct Same
    {
        static constexpr bool value = std::is_same_v<T,U1> && Same<T,Us...>::value;
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

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        template <class...> using Void_t = void;

        struct Nonesuch {
            Nonesuch() = delete;
            ~Nonesuch() = delete;
            Nonesuch(Nonesuch const&) = delete;
            Nonesuch(Nonesuch &&) = delete;
            void operator=(Nonesuch const&) = delete;
            void operator=(Nonesuch &&) = delete;
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
    /// \endcond

    /// \ingroup amrex_utilities
    //! True if the expression Op<Args...> is valid (SFINAE detection idiom).
    template <template <class...> class Op, class... Args>
    using IsDetected = typename detail::Detector<detail::Nonesuch, void, Op, Args...>::value_t;

    /// \ingroup amrex_utilities
    //! The type Op<Args...> if valid, otherwise detail::Nonesuch.
    template <template <class...> class Op, class... Args>
    using Detected_t = typename detail::Detector<detail::Nonesuch, void, Op, Args...>::type;

    /// \ingroup amrex_utilities
    //! The type Op<Args...> if valid, otherwise Default.
    template <class Default, template <class...> class Op, class... Args>
    using DetectedOr = typename detail::Detector<Default, void, Op, Args...>::type;

    /// \ingroup amrex_utilities
    //! True if Op<Args...> is valid and equals Expected.
    template <class Expected, template <typename...> class Op, class... Args>
    using IsDetectedExact = std::is_same<Expected, Detected_t<Op, Args...>>;

    ////////////////////////////////////////////////////////////////////////////////
    //                                                           [traits.IsCallable]

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        template <typename T, typename... Args>
        using call_result_t = decltype(std::declval<T>()(std::declval<Args>()...));
    }
    /// \endcond

    //! \ingroup amrex_utilities
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

    //! \ingroup amrex_utilities
    //! \brief Test if a given type T is callable with arguments of type Args... and returns R.
    //!
    //! \see IsCallable
    template <typename R, typename T, typename... Args>
    struct IsCallableR : IsDetectedExact<R, detail::call_result_t, T, Args...> {};

    ////////////////////////////////////////////////////////////////////////////////
    //                                                          [traits.Conjunction]
    //                                                          [traits.Disjunction]
    //                                                          [traits.Negation]

#if defined(__cpp_lib_logical_traits)
    /// \ingroup amrex_utilities
    //! Logical AND of compile-time boolean predicates (equivalent to std::conjunction).
    template <typename... Args> using Conjunction = std::conjunction<Args...>;
    /// \ingroup amrex_utilities
    //! Logical OR of compile-time boolean predicates (equivalent to std::disjunction).
    template <typename... Args> using Disjunction = std::disjunction<Args...>;
    /// \ingroup amrex_utilities
    //! Logical NOT of a compile-time boolean predicate (equivalent to std::negation).
    template <typename... Args> using Negation = std::negation<Args...>;
/// \cond DOXYGEN_IGNORE
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
/// \endcond

    ////////////////////////////////////////////////////////////////////////////////
    //                                                        [traits.IsConvertible]

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        template<typename T>
        inline constexpr bool is_convertible(T) {return true;}

        template <typename T, typename U, typename Enable = void>
        struct IsConvertibleImp : std::false_type {};

        template <typename T, typename U>
        struct IsConvertibleImp<T, U, std::enable_if_t<is_convertible<T>(U{})>> : std::true_type {};
    }
    /// \endcond

    //! \ingroup amrex_utilities
    //! \brief Test if all the types Args... are automatically convertible to type T.
    template <typename T, typename... Args>
    struct IsConvertible {
        static constexpr bool value = (... && detail::IsConvertibleImp<T, Args>::value);
    };

    /// \ingroup amrex_utilities
    //! True if all types Args... are automatically convertible to T.
    template <typename T, typename... Args>
    inline constexpr bool IsConvertible_v = IsConvertible<T, Args...>::value;

    /// \cond DOXYGEN_IGNORE
    namespace detail {
        template <typename T, typename U, typename Enable = void>
        struct IsNarrowingConversionImp : std::true_type {};

        template <typename T, typename U>
        struct IsNarrowingConversionImp<T, U, std::void_t<decltype(U{std::declval<T>()})>> : std::false_type {};
    }
    /// \endcond

    /**
     * \class amrex::IsNarrowingConversion
     * \ingroup amrex_utilities
     * \brief True if converting From to To is a narrowing conversion.
     */
    template <typename From, typename To>
    struct IsNarrowingConversion : detail::IsNarrowingConversionImp<From, To> {};

    /// \ingroup amrex_utilities
    //! True if converting From to To is a narrowing conversion.
    template <typename From, typename To>
    inline constexpr bool IsNarrowingConversion_v = IsNarrowingConversion<From, To>::value;

    /// \ingroup amrex_utilities
    //! True if converting From to To is NOT a narrowing conversion.
    template <typename From, typename To>
    inline constexpr bool IsNonNarrowingConversion_v = !IsNarrowingConversion<From, To>::value;

    // Move this down, because doxygen can not parse anything below IsStoreAtomic
    /**
     * \class amrex::IsStoreAtomic
     * \ingroup amrex_utilities
     * \brief True if T is an arithmetic type that can be stored atomically (8 bytes or fewer).
     */
    template <class T, class Enable = void>
    struct IsStoreAtomic : std::false_type {};
    //
    template <class T>
    struct IsStoreAtomic<T, std::enable_if_t <
                                std::is_arithmetic_v<T>
                                && sizeof(T) <= 8 > >
        : std::true_type {};

    /**
     * \class amrex::IsStdVector
     * \ingroup amrex_utilities
     * \brief True if T is derived from std::vector<T::value_type, T::allocator_type>.
     */
    template <class T, class Enable = void>
    struct IsStdVector : std::false_type {};
    //
    template <class T>
    struct IsStdVector<T, std::enable_if_t<std::is_base_of_v<std::vector<typename T::value_type,
                                                                       typename T::allocator_type>,
                                                           T>> >
                       : std::true_type {};

}

#endif
