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
concept BaseFabType = amrex::IsBaseFab<T>::value;

/**
 * \brief Checks if a type is an amrex::FabArray (like MultiFab).
 */
template <typename T>
concept FabArrayType = amrex::IsFabArray<T>::value;  
} // namespace amrex
#endif
