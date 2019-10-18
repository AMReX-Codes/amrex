#ifndef WARPX_amrex_qed_wrapper_commons_h_
#define WARPX_amrex_qed_wrapper_commons_h_

//Common definitions for the QED library wrappers and table builders

#include <AMReX_AmrCore.H>
#include <AMReX_GpuQualifiers.H>

//Sets the decorator for GPU
#define PXRMP_GPU AMREX_GPU_DEVICE
//Sets SI units in the library
#define PXRMP_WITH_SI_UNITS

namespace QedUtils{
    //An empty data type
    struct DummyStruct{};
};

//Control structures used by the engines
struct PicsarBreitWheelerCtrl
{
              _REAL chi_phot_min =
            static_cast<_REAL>(__breit_wheeler_min_chi_phot);

          _REAL chi_phot_tdndt_min =
            static_cast<_REAL>(__breit_wheeler_min_tdndt_chi_phot);
          _REAL chi_phot_tdndt_max =
            static_cast<_REAL>(__breit_wheeler_max_tdndt_chi_phot);
          size_t chi_phot_tdndt_how_many =
            __breit_wheeler_how_many_tdndt_chi_phot;

          _REAL chi_phot_tpair_min =
            static_cast<_REAL>(__breit_wheeler_min_tpair_chi_phot);
          _REAL chi_phot_tpair_max =
            static_cast<_REAL>(__breit_wheeler_max_tpair_chi_phot);
          size_t chi_phot_tpair_how_many =
            __breit_wheeler_how_many_tpair_chi_phot;
          size_t chi_frac_tpair_how_many =
            __breit_wheeler_chi_frac_tpair_how_many;
};

#endif //WARPX_amrex_qed_wrapper_commons_h_
