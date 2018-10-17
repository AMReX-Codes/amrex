#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_

#include <type_traits>

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

    template <class A, class Enable = void> struct IsTriviallyCopyable : std::false_type {};

    // workaround for the fact that GCC < 5 doesn't have std::is_trivially_copyable
    template <class D>
    struct IsTriviallyCopyable<D, typename std::enable_if<
#if __GNUG__ && __GNUC__ < 5    
                           __has_trivial_copy(D) > >
#else
                           std::is_trivially_copyable<D>::value> >
#endif
        : std::true_type {};

}

#endif
