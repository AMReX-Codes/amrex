#ifndef AMREX_FAB_DATA_TYPE_H_
#define AMREX_FAB_DATA_TYPE_H_
#include <AMReX_Config.H>

#include <AMReX_TypeTraits.H>

namespace amrex {

template <typename T> struct FabDataType {};
//
template <typename T>
    requires(IsMultiFabLike_v<T>)
struct FabDataType<T>
{
    using   fab_type = typename T::fab_type;
    using value_type = typename T::value_type;
};

template <typename T>
    requires(IsMultiFabLike_v<typename T::value_type>)
struct FabDataType<T>
{
    using   fab_type = typename T::value_type::fab_type;
    using value_type = typename T::value_type::value_type;
};

}

#endif
