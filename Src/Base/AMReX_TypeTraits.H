#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_

#include <type_traits>

namespace amrex
{
    // BaseFab
    template <class T> class BaseFab;
    //
    template <class A, class Enable = void> struct IsBaseFab : std::false_type {};
    //
    template <class D>
    struct IsBaseFab<D, typename std::enable_if<std::is_base_of<BaseFab<typename D::value_type>,D>::value>::type> : std::true_type {};

}

#endif
