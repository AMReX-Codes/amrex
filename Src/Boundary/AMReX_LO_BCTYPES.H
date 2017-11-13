#ifndef _LO_BCTYPES_H_
#define _LO_BCTYPES_H_

#define LO_DIRICHLET 101
#define LO_NEUMANN 102
#define LO_REFLECT_ODD 103
#define LO_MARSHAK 104
#define LO_SANCHEZ_POMRANING 105

#ifdef __cplusplus
namespace amrex {
    enum struct LinOpBCType : int {
        interior         = 0,
        Dirichlet        = LO_DIRICHLET,
        Neumann          = LO_NEUMANN,
        reflect_odd      = LO_REFLECT_ODD,
        Marshak          = LO_MARSHAK,
        SanchezPomraning = LO_SANCHEZ_POMRANING,
        Periodic,
        bogus    
    };
}
#endif

#endif

