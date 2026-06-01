#ifndef AMREX_CONCEPTS_H_
#define AMREX_CONCEPTS_H_

#include <AMReX_Config.H>
#include <AMReX_TypeTraits.H>
#include <concepts>

namespace amrex {
/**
 * \brief Checks if a type is derived from amrex::BaseFab.
 */
template <typename T>
concept BaseFabType = amrex::IsBaseFab_v<T>;

/**
 * \brief Checks if a type is an amrex::FabArray.
 */
template <typename T>
concept FabArrayType = amrex::IsFabArray_v<T>;

/**
 * \brief Checks if a type is amrex::MultiFab-like (i.e., FabArray<FAB>,
 * where FAB is a BaseFabType).
 */
template <typename T>
concept MultiFabLike = amrex::IsMultiFabLike_v<T>;

} // namespace amrex
#endif
