#ifndef AMREX_EB_MULTIFAB_UTIL_3D_C_H_
#define AMREX_EB_MULTIFAB_UTIL_3D_C_H_

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v;
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void eb_set_covered_nodes (int i, int j, int k, int n, int icomp, Array4<Real> const& d,
                           Array4<EBCellFlag const> const& f, Real const * AMREX_RESTRICT v)
{
    if (f(i-1,j-1,k-1).isCovered() and f(i  ,j-1,k-1).isCovered() and
        f(i-1,j  ,k-1).isCovered() and f(i  ,j  ,k-1).isCovered() and
        f(i-1,j-1,k  ).isCovered() and f(i  ,j-1,k  ).isCovered() and
        f(i-1,j  ,k  ).isCovered() and f(i  ,j  ,k  ).isCovered())
    {
        d(i,j,k,n+icomp) = v[n];
    }
}

}

#endif
