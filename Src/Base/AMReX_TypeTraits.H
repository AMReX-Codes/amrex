#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_

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
    // CUDA supports atomicAdd for the following types
    template <> struct HasAtomicAdd<int> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned int> : std::true_type {};
    template <> struct HasAtomicAdd<unsigned long long> : std::true_type {};
    template <> struct HasAtomicAdd<float> : std::true_type {};
    template <> struct HasAtomicAdd<double> : std::true_type {};

#ifdef AMREX_PARTICLES
    template <bool is_const, int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    class ParIterBase;

    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
    class ParIter;

    template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
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
        static constexpr bool value = std::is_same<T,U1>::value and Same<T,Us...>::value;
    };
}

#endif
