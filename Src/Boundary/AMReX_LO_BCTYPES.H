#ifndef AMREX_LO_BCTYPES_H_
#define AMREX_LO_BCTYPES_H_

#define LO_DIRICHLET 101
#define LO_NEUMANN 102
#define LO_REFLECT_ODD 103
#define LO_MARSHAK 104
#define LO_SANCHEZ_POMRANING 105

#define AMREX_LO_DIRICHLET 101
#define AMREX_LO_NEUMANN 102
#define AMREX_LO_REFLECT_ODD 103
#define AMREX_LO_MARSHAK 104
#define AMREX_LO_SANCHEZ_POMRANING 105
#define AMREX_LO_INFLOW   106
#define AMREX_LO_PERIODIC 200
#define AMREX_LO_BOGUS    1729

#ifdef __cplusplus
namespace amrex {
    enum struct LinOpBCType : int {
        interior         = 0,
        Dirichlet        = LO_DIRICHLET,
        Neumann          = LO_NEUMANN,
        reflect_odd      = LO_REFLECT_ODD,
        Marshak          = LO_MARSHAK,
        SanchezPomraning = LO_SANCHEZ_POMRANING,
        inflow           = AMREX_LO_INFLOW,
        Periodic         = AMREX_LO_PERIODIC,
        bogus            = AMREX_LO_BOGUS
    };
}
#endif

#endif

