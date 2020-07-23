#include <AMReX_FilND_C.H>

namespace amrex {

void fab_filnd (Box const& bx, Array4<Real> const& qn, int ncomp,
                Box const& domain, Real const* /*dx*/, Real const* /*xlo*/,
                BCRec const* bcn)
{
    const auto lo = amrex::lbound(bx);
    const auto hi = amrex::ubound(bx);
    const auto domlo = amrex::lbound(domain);
    const auto domhi = amrex::ubound(domain);

    const int ilo = domlo.x;
    const int ihi = domhi.x;

#if AMREX_SPACEDIM >= 2
    const int jlo = domlo.y;
    const int jhi = domhi.y;
#endif

#if AMREX_SPACEDIM == 3
    const int klo = domlo.z;
    const int khi = domhi.z;
#endif

    for (int n = 0; n < ncomp; ++n)
    {
        Array4<Real> q(qn,n);
        BCRec const& bc = bcn[n];

        if (lo.x < ilo && (bc.lo(0) != BCType::int_dir)) {
           const int imin = lo.x;
           const int imax = ilo-1;
           for (int k = lo.z; k <= hi.z; ++k) {
           for (int j = lo.y; j <= hi.y; ++j) {
           for (int i = imin; i <= imax; ++i) {
               q(i,j,k) = q(ilo,j,k);
           }}}
        }

        if (hi.x > ihi && (bc.hi(0) != BCType::int_dir)) {
            const int imin = ihi+1;
            const int imax = hi.x;
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = imin; i <= imax; ++i) {
                q(i,j,k) = q(ihi,j,k);
            }}}
        }

#if AMREX_SPACEDIM >= 2

        if (lo.y < jlo && (bc.lo(1) != BCType::int_dir)) {
            const int jmin = lo.y;
            const int jmax = jlo-1;
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = jmin; j <= jmax; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                q(i,j,k) = q(i,jlo,k);
            }}}
        }

        if (hi.y > jhi && (bc.hi(1) != BCType::int_dir)) {
            const int jmin = jhi+1;
            const int jmax = hi.y;
            for (int k = lo.z; k <= hi.z; ++k) {
            for (int j = jmin; j <= jmax; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                q(i,j,k) = q(i,jhi,k);
            }}}
        }
#endif

#if AMREX_SPACEDIM == 3

        if (lo.z < klo && (bc.lo(2) != BCType::int_dir)) {
            const int kmin = lo.z;
            const int kmax = klo-1;
            for (int k = kmin; k <= kmax; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                q(i,j,k) = q(i,j,klo);
            }}}
        }

        if (hi.z > khi && (bc.hi(2) != BCType::int_dir)) {
            const int kmin = khi+1;
            const int kmax = hi.z;
            for (int k = kmin; k <= kmax; ++k) {
            for (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                q(i,j,k) = q(i,j,khi);
            }}}
        }
#endif
    }
}

}
