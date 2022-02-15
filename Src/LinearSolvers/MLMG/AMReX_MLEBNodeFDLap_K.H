#ifndef AMREX_MLEBNODEFDLAP_K_H_
#define AMREX_MLEBNODEFDLAP_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLEBNodeFDLaplacian.H>
#include <AMReX_LO_BCTYPES.H>

#if (AMREX_SPACEDIM == 1)
#include <AMReX_MLEBNodeFDLap_1D_K.H>
#elif (AMREX_SPACEDIM == 2)
#include <AMReX_MLEBNodeFDLap_2D_K.H>
#else
#include <AMReX_MLEBNodeFDLap_3D_K.H>
#endif

namespace amrex {

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x_doit (int i, int j, int k, Array4<Real> const& px,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecx, F&& phieb, Real dxi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i+1,j,k) >= 0) {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i+1,j,k) < 0) {
        px(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        px(i,j,k) = dxi * (p(i+1,j,k) - phieb(i,j,k)) / (Real(1.0) - Real(2.0) * ecx(i,j,k));
    } else { //
        px(i,j,k) = dxi * (phieb(i+1,j,k) - p(i,j,k)) / (Real(1.0) + Real(2.0) * ecx(i,j,k));
    }
}

template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y_doit (int i, int j, int k, Array4<Real> const& py,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecy, F&& phieb, Real dyi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j+1,k) >= 0) {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j+1,k) < 0) {
        py(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        py(i,j,k) = dyi * (p(i,j+1,k) - phieb(i,j,k)) / (Real(1.0) - Real(2.0) * ecy(i,j,k));
    } else { //
        py(i,j,k) = dyi * (phieb(i,j+1,k) - p(i,j,k)) / (Real(1.0) + Real(2.0) * ecy(i,j,k));
    }
}

#if (AMREX_SPACEDIM > 2)
template <typename F>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z_doit (int i, int j, int k, Array4<Real> const& pz,
                              Array4<Real const> const& p, Array4<int const> const& dmsk,
                              Array4<Real const> const& ecz, F&& phieb, Real dzi)
{
    if (dmsk(i,j,k) >= 0 && dmsk(i,j,k+1) >= 0) {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    } else if (dmsk(i,j,k) < 0 && dmsk(i,j,k+1) < 0) {
        pz(i,j,k) = Real(0.0);
    } else if (dmsk(i,j,k) < 0) {
        pz(i,j,k) = dzi * (p(i,j,k+1) - phieb(i,j,k)) / (Real(1.0) - Real(2.0) * ecz(i,j,k));
    } else { //
        pz(i,j,k) = dzi * (phieb(i,j,k+1) - p(i,j,k)) / (Real(1.0) + Real(2.0) * ecz(i,j,k));
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                         Real phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ecx,
                                [=] (int, int, int) -> Real { return phieb; },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecx,
                         Array4<Real const> const& phieb, Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_x_doit(i,j,k, px, p, dmsk, ecx,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dxi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecy,
                         Real phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ecy,
                                [=] (int, int, int) -> Real { return phieb; },
                                dyi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecy,
                         Array4<Real const> const& phieb, Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_y_doit(i,j,k, py, p, dmsk, ecy,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dyi);
    });
}

#if (AMREX_SPACEDIM > 2)

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecz,
                         Real phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ecz,
                                [=] (int, int, int) -> Real { return phieb; },
                                dzi);
    });
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Array4<int const> const& dmsk, Array4<Real const> const& ecz,
                         Array4<Real const> const& phieb, Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        mlebndfdlap_grad_z_doit(i,j,k, pz, p, dmsk, ecz,
                                [=] (int i1, int i2, int i3) -> Real { return phieb(i1,i2,i3); },
                                dzi);
    });
}

#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_x (Box const& b, Array4<Real> const& px, Array4<Real const> const& p,
                         Real dxi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        px(i,j,k) = dxi * (p(i+1,j,k) - p(i,j,k));
    });
}

#if (AMREX_SPACEDIM > 1)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_y (Box const& b, Array4<Real> const& py, Array4<Real const> const& p,
                         Real dyi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        py(i,j,k) = dyi * (p(i,j+1,k) - p(i,j,k));
    });
}
#endif

#if (AMREX_SPACEDIM > 2)
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mlebndfdlap_grad_z (Box const& b, Array4<Real> const& pz, Array4<Real const> const& p,
                         Real dzi)
{
    AMREX_LOOP_3D(b, i, j, k,
    {
        pz(i,j,k) = dzi * (p(i,j,k+1) - p(i,j,k));
    });
}
#endif

}

#endif
