// Do not include this header anywhere other than AMReX_Loop.H.
// The purpose of this file is to avoid clang-tidy.

#define AMREX_LOOP_3D(bx, i, j, k, block) \
{ \
    const auto amrex_i_lo = amrex::lbound(bx); \
    const auto amrex_i_hi = amrex::ubound(bx); \
    for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
    for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
    for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
        AMREX_D_PICK(amrex::ignore_unused(j,k),amrex::ignore_unused(k),(void)0); \
        block \
    }}} \
}

#define AMREX_LOOP_4D(bx, ncomp, i, j, k, n, block) \
{ \
    const auto amrex_i_lo = amrex::lbound(bx); \
    const auto amrex_i_hi = amrex::ubound(bx); \
    for (int n = 0; n < ncomp; ++n) { \
    for (int k = amrex_i_lo.z; k <= amrex_i_hi.z; ++k) { \
    for (int j = amrex_i_lo.y; j <= amrex_i_hi.y; ++j) { \
    for (int i = amrex_i_lo.x; i <= amrex_i_hi.x; ++i) { \
        AMREX_D_PICK(amrex::ignore_unused(j,k),amrex::ignore_unused(k),(void)0); \
        block \
    }}}} \
}
