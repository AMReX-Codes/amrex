#ifndef AMREX_CONCEPTS_H_
#define AMREX_CONCEPTS_H_

#include <AMReX_Config.H>
#include <AMReX_TypeTraits.H>
#include <concepts>

namespace amrex {
  
  /**
   * \brief Checks if a type is derived from amrex::BaseFab.
   * Matches your requirement: requires IsBaseFabType<MFFAB>
   */
  template <typename T>
  concept IsBaseFabType = amrex::IsBaseFab<T>::value;
  
  /**
   * \brief Checks if a type is an amrex::FabArray (like MultiFab).
   * Matches your requirement: requires IsFabArrayType<MFFAB>
   */
  template <typename T>
  concept IsFabArrayType = amrex::IsFabArray<T>::value;

  template <typename T, typename Interp>
  concept IsFabArrayTypeWithSpecificInterp = amrex::IsFabArray<T>::value && (!std::is_same_v<Interp, amrex::MFInterpolater>);
  
  
} // namespace amrex
#endif
