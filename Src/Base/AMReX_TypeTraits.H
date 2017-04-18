#ifndef BL_TYPETRAITS_H_
#define BL_TYPETRAITS_H_

#include <type_traits>

namespace amrex
{
    template <class A, class Enable = void> struct IsBaseFab : std::false_type {};
}

#endif
