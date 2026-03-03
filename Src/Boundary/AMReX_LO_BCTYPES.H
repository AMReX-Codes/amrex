#ifndef AMREX_LO_BCTYPES_H_
#define AMREX_LO_BCTYPES_H_
#include <AMReX_Config.H>

#define AMREX_LO_DIRICHLET 101
#define AMREX_LO_NEUMANN 102
#define AMREX_LO_REFLECT_ODD 103
#define AMREX_LO_MARSHAK 104
#define AMREX_LO_SANCHEZ_POMRANING 105
#define AMREX_LO_INFLOW   106
#define AMREX_LO_INHOMOG_NEUMANN 107
#define AMREX_LO_ROBIN 108
#define AMREX_LO_SYMMETRY 109
#define AMREX_LO_PERIODIC 200
#define AMREX_LO_BOGUS    1729

/**
 * \file AMReX_LO_BCTYPES.H
 * \brief Constants describing linear-operator boundary condition types.
 */

#ifdef __cplusplus

#include <iosfwd>

namespace amrex {
    enum struct LinOpBCType : int {
        interior         = 0,
        Dirichlet        = AMREX_LO_DIRICHLET,
        Neumann          = AMREX_LO_NEUMANN,
        reflect_odd      = AMREX_LO_REFLECT_ODD,
        Marshak          = AMREX_LO_MARSHAK,
        SanchezPomraning = AMREX_LO_SANCHEZ_POMRANING,
        inflow           = AMREX_LO_INFLOW,
        inhomogNeumann   = AMREX_LO_INHOMOG_NEUMANN,
        Robin            = AMREX_LO_ROBIN,
        symmetry         = AMREX_LO_SYMMETRY,
        Periodic         = AMREX_LO_PERIODIC,
        bogus            = AMREX_LO_BOGUS
        // If we add a new type, we should update AMReX_lo_bctypes_mode.F90.
    };

    /**
     * \brief Stream boundary-condition type names for debugging.
     *
     * \param os Output stream to receive the text label.
     * \param t  Boundary-condition type to print.
     * \return The output stream for chaining.
     */
    std::ostream& operator<< (std::ostream& os, const amrex::LinOpBCType& t);
}

#endif

#endif
