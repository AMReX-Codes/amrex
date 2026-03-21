#ifndef AMREX_FAB_DATA_TYPE_H_
#define AMREX_FAB_DATA_TYPE_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_FabDataType.H
 * \brief Type trait that exposes the FAB and value types of MultiFab-like containers.
 */

#include <AMReX_Concepts.H>
#include <AMReX_TypeTraits.H>

namespace amrex {

//! Type trait specialized for MultiFab-like types and their containers.
template <typename T> struct FabDataType {};
//
//! Extract the FAB and scalar types directly from a MultiFab-like type.
template <MultiFabLike T>
struct FabDataType<T>
{
    //! FAB type stored by the MultiFab-like type.
    using   fab_type = typename T::fab_type;
    //! Scalar value type held by the FAB.
    using value_type = typename T::value_type;
};

//! Extract the FAB and scalar types when the container element is MultiFab-like.
template <typename T>
requires (MultiFabLike<typename T::value_type>)
struct FabDataType<T>
{
    //! FAB type stored by the container's MultiFab-like element.
    using   fab_type = typename T::value_type::fab_type;
    //! Scalar value type held by that FAB.
    using value_type = typename T::value_type::value_type;
};

}

#endif
