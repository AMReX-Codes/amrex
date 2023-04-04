#ifndef AMREX_MLTENSOR_3D_K_H_
#define AMREX_MLTENSOR_3D_K_H_
#include <AMReX_Config.H>

#include <AMReX_MLLinOp_K.H>

namespace amrex {

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xlo_ylo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxlo,
                                  Array4<int const> const& mylo,
                                  Array4<Real const> const& bcvalxlo,
                                  Array4<Real const> const& bcvalylo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xlo_domain, bool ylo_domain) noexcept
{
    if (mxlo(i,j,k) != BndryData::covered && (!xlo_domain || !ylo_domain)) {
        bool x_interior = mylo(i+1,j  ,k) == BndryData::covered;
        bool x_exterior = mylo(i+1,j  ,k) == BndryData::not_covered;
        bool y_interior = mxlo(i  ,j+1,k) == BndryData::covered;
        bool y_exterior = mxlo(i  ,j+1,k) == BndryData::not_covered;
        if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (y_interior || ylo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xhi_ylo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxhi,
                                  Array4<int const> const& mylo,
                                  Array4<Real const> const& bcvalxhi,
                                  Array4<Real const> const& bcvalylo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xhi_domain, bool ylo_domain) noexcept
{
    if (mxhi(i,j,k) != BndryData::covered && (!xhi_domain || !ylo_domain)) {
        bool x_interior = mylo(i-1,j  ,k) == BndryData::covered;
        bool x_exterior = mylo(i-1,j  ,k) == BndryData::not_covered;
        bool y_interior = mxhi(i  ,j+1,k) == BndryData::covered;
        bool y_exterior = mxhi(i  ,j+1,k) == BndryData::not_covered;
        if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (y_interior || ylo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xlo_yhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxlo,
                                  Array4<int const> const& myhi,
                                  Array4<Real const> const& bcvalxlo,
                                  Array4<Real const> const& bcvalyhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xlo_domain, bool yhi_domain) noexcept
{
    if (mxlo(i,j,k) != BndryData::covered && (!xlo_domain || !yhi_domain)) {
        bool x_interior = myhi(i+1,j  ,k) == BndryData::covered;
        bool x_exterior = myhi(i+1,j  ,k) == BndryData::not_covered;
        bool y_interior = mxlo(i  ,j-1,k) == BndryData::covered;
        bool y_exterior = mxlo(i  ,j-1,k) == BndryData::not_covered;
        if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (y_interior || yhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xhi_yhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxhi,
                                  Array4<int const> const& myhi,
                                  Array4<Real const> const& bcvalxhi,
                                  Array4<Real const> const& bcvalyhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xhi_domain, bool yhi_domain) noexcept
{
    if (mxhi(i,j,k) != BndryData::covered && (!xhi_domain || !yhi_domain)) {
        bool x_interior = myhi(i-1,j  ,k) == BndryData::covered;
        bool x_exterior = myhi(i-1,j  ,k) == BndryData::not_covered;
        bool y_interior = mxhi(i  ,j-1,k) == BndryData::covered;
        bool y_exterior = mxhi(i  ,j-1,k) == BndryData::not_covered;
        if ((x_interior && y_interior) || (x_exterior && y_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (y_interior || yhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xlo_zlo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxlo,
                                  Array4<int const> const& mzlo,
                                  Array4<Real const> const& bcvalxlo,
                                  Array4<Real const> const& bcvalzlo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xlo_domain, bool zlo_domain) noexcept
{
    if (mxlo(i,j,k) != BndryData::covered && (!xlo_domain || !zlo_domain)) {
        bool x_interior = mzlo(i+1,j,k  ) == BndryData::covered;
        bool x_exterior = mzlo(i+1,j,k  ) == BndryData::not_covered;
        bool z_interior = mxlo(i  ,j,k+1) == BndryData::covered;
        bool z_exterior = mxlo(i  ,j,k+1) == BndryData::not_covered;
        if ((x_interior && z_interior) || (x_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (z_interior || zlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xhi_zlo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxhi,
                                  Array4<int const> const& mzlo,
                                  Array4<Real const> const& bcvalxhi,
                                  Array4<Real const> const& bcvalzlo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xhi_domain, bool zlo_domain) noexcept
{
    if (mxhi(i,j,k) != BndryData::covered && (!xhi_domain || !zlo_domain)) {
        bool x_interior = mzlo(i-1,j,k  ) == BndryData::covered;
        bool x_exterior = mzlo(i-1,j,k  ) == BndryData::not_covered;
        bool z_interior = mxhi(i  ,j,k+1) == BndryData::covered;
        bool z_exterior = mxhi(i  ,j,k+1) == BndryData::not_covered;
        if ((x_interior && z_interior) || (x_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (z_interior || zlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xlo_zhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxlo,
                                  Array4<int const> const& mzhi,
                                  Array4<Real const> const& bcvalxlo,
                                  Array4<Real const> const& bcvalzhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xlo_domain, bool zhi_domain) noexcept
{
    if (mxlo(i,j,k) != BndryData::covered && (!xlo_domain || !zhi_domain)) {
        bool x_interior = mzhi(i+1,j,k  ) == BndryData::covered;
        bool x_exterior = mzhi(i+1,j,k  ) == BndryData::not_covered;
        bool z_interior = mxlo(i  ,j,k-1) == BndryData::covered;
        bool z_exterior = mxlo(i  ,j,k-1) == BndryData::not_covered;
        if ((x_interior && z_interior) || (x_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                   bct(Orientation::xlo(), icomp),
                                   bcl(Orientation::xlo(), icomp),
                                   bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (z_interior || zhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_xhi_zhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mxhi,
                                  Array4<int const> const& mzhi,
                                  Array4<Real const> const& bcvalxhi,
                                  Array4<Real const> const& bcvalzhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool xhi_domain, bool zhi_domain) noexcept
{
    if (mxhi(i,j,k) != BndryData::covered && (!xhi_domain || !zhi_domain)) {
        bool x_interior = mzhi(i-1,j,k  ) == BndryData::covered;
        bool x_exterior = mzhi(i-1,j,k  ) == BndryData::not_covered;
        bool z_interior = mxhi(i  ,j,k-1) == BndryData::covered;
        bool z_exterior = mxhi(i  ,j,k-1) == BndryData::not_covered;
        if ((x_interior && z_interior) || (x_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (x_interior || xhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                   bct(Orientation::xhi(), icomp),
                                   bcl(Orientation::xhi(), icomp),
                                   bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
            }
        } else if (z_interior || zhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_ylo_zlo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mylo,
                                  Array4<int const> const& mzlo,
                                  Array4<Real const> const& bcvalylo,
                                  Array4<Real const> const& bcvalzlo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool ylo_domain, bool zlo_domain) noexcept
{
    if (mylo(i,j,k) != BndryData::covered && (!ylo_domain || !zlo_domain)) {
        bool y_interior = mzlo(i,j+1,k  ) == BndryData::covered;
        bool y_exterior = mzlo(i,j+1,k  ) == BndryData::not_covered;
        bool z_interior = mylo(i,j  ,k+1) == BndryData::covered;
        bool z_exterior = mylo(i,j  ,k+1) == BndryData::not_covered;
        if ((y_interior && z_interior) || (y_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (y_interior || ylo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
            }
        } else if (z_interior || zlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_yhi_zlo (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& myhi,
                                  Array4<int const> const& mzlo,
                                  Array4<Real const> const& bcvalyhi,
                                  Array4<Real const> const& bcvalzlo,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool yhi_domain, bool zlo_domain) noexcept
{
    if (myhi(i,j,k) != BndryData::covered && (!yhi_domain || !zlo_domain)) {
        bool y_interior = mzlo(i,j-1,k  ) == BndryData::covered;
        bool y_exterior = mzlo(i,j-1,k  ) == BndryData::not_covered;
        bool z_interior = myhi(i,j  ,k+1) == BndryData::covered;
        bool z_exterior = myhi(i,j  ,k+1) == BndryData::not_covered;
        if ((y_interior && z_interior) || (y_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (y_interior || yhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
            }
        } else if (z_interior || zlo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                   bct(Orientation::zlo(), icomp),
                                   bcl(Orientation::zlo(), icomp),
                                   bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_ylo_zhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& mylo,
                                  Array4<int const> const& mzhi,
                                  Array4<Real const> const& bcvalylo,
                                  Array4<Real const> const& bcvalzhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool ylo_domain, bool zhi_domain) noexcept
{
    if (mylo(i,j,k) != BndryData::covered && (!ylo_domain || !zhi_domain)) {
        bool y_interior = mzhi(i,j+1,k  ) == BndryData::covered;
        bool y_exterior = mzhi(i,j+1,k  ) == BndryData::not_covered;
        bool z_interior = mylo(i,j  ,k-1) == BndryData::covered;
        bool z_exterior = mylo(i,j  ,k-1) == BndryData::not_covered;
        if ((y_interior && z_interior) || (y_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (y_interior || ylo_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                   bct(Orientation::ylo(), icomp),
                                   bcl(Orientation::ylo(), icomp),
                                   bcvalylo, maxorder, dxinv[1], inhomog, icomp);
            }
        } else if (z_interior || zhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges_yhi_zhi (int const i, int const j, int const k, Dim3 const& blen,
                                  Array4<Real> const& vel,
                                  Array4<int const> const& myhi,
                                  Array4<int const> const& mzhi,
                                  Array4<Real const> const& bcvalyhi,
                                  Array4<Real const> const& bcvalzhi,
                                  Array2D<BoundCond,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bct,
                                  Array2D<Real,
                                          0,2*AMREX_SPACEDIM,
                                          0,AMREX_SPACEDIM> const& bcl,
                                  int inhomog, int maxorder,
                                  GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                                  bool yhi_domain, bool zhi_domain) noexcept
{
    if (myhi(i,j,k) != BndryData::covered && (!yhi_domain || !zhi_domain)) {
        bool y_interior = mzhi(i,j-1,k  ) == BndryData::covered;
        bool y_exterior = mzhi(i,j-1,k  ) == BndryData::not_covered;
        bool z_interior = myhi(i,j  ,k-1) == BndryData::covered;
        bool z_exterior = myhi(i,j  ,k-1) == BndryData::not_covered;
        if ((y_interior && z_interior) || (y_exterior && z_exterior)) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                Real tmp = vel(i,j,k,icomp);
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
            }
        } else if (y_interior || yhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                   bct(Orientation::yhi(), icomp),
                                   bcl(Orientation::yhi(), icomp),
                                   bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
            }
        } else if (z_interior || zhi_domain) {
            for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
                mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                   bct(Orientation::zhi(), icomp),
                                   bcl(Orientation::zhi(), icomp),
                                   bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
            }
        }
    }
}

#ifdef AMREX_USE_EB
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_corners (int icorner, Box const& vbox, // vbox: the valid box
                            Array4<Real> const& vel,
                            Array4<int const> const& mxlo,
                            Array4<int const> const& mylo,
                            Array4<int const> const& mzlo,
                            Array4<int const> const& mxhi,
                            Array4<int const> const& myhi,
                            Array4<int const> const& mzhi,
                            Array4<Real const> const& bcvalxlo,
                            Array4<Real const> const& bcvalylo,
                            Array4<Real const> const& bcvalzlo,
                            Array4<Real const> const& bcvalxhi,
                            Array4<Real const> const& bcvalyhi,
                            Array4<Real const> const& bcvalzhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Array2D<Real,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bcl,
                            int inhomog, int maxorder,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    bool xlo_domain = (vlo.x == dlo.x);
    bool ylo_domain = (vlo.y == dlo.y);
    bool zlo_domain = (vlo.z == dlo.z);
    bool xhi_domain = (vhi.x == dhi.x);
    bool yhi_domain = (vhi.y == dhi.y);
    bool zhi_domain = (vhi.z == dhi.z);

    for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
        switch (icorner) {
        case 0: {
            // xlo & ylo & zlo
            int i = vlo.x-1;
            int j = vlo.y-1;
            int k = vlo.z-1;
            if (mxlo(i,j,k) != BndryData::covered &&
                (!xlo_domain || !ylo_domain || !zlo_domain)) {
                bool x_interior = mylo(i+1,j  ,k  ) == BndryData::covered;
                bool x_exterior = mylo(i+1,j  ,k  ) == BndryData::not_covered;
                bool y_interior = mxlo(i  ,j+1,k  ) == BndryData::covered;
                bool y_exterior = mxlo(i  ,j+1,k  ) == BndryData::not_covered;
                bool z_interior = mxlo(i  ,j  ,k+1) == BndryData::covered;
                bool z_exterior = mxlo(i  ,j  ,k+1) == BndryData::not_covered;
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 1: {
            // xhi & ylo & zlo
            int i = vhi.x+1;
            int j = vlo.y-1;
            int k = vlo.z-1;
            bool x_interior = mylo(i-1,j  ,k  ) == BndryData::covered;
            bool x_exterior = mylo(i-1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j+1,k  ) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j+1,k  ) == BndryData::not_covered;
            bool z_interior = mxhi(i  ,j  ,k+1) == BndryData::covered;
            bool z_exterior = mxhi(i  ,j  ,k+1) == BndryData::not_covered;
            if (mxhi(i,j,k) != BndryData::covered &&
                (!xhi_domain || !ylo_domain || !zlo_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 2: {
            // xlo & yhi & zlo
            int i = vlo.x-1;
            int j = vhi.y+1;
            int k = vlo.z-1;
            bool x_interior = myhi(i+1,j  ,k  ) == BndryData::covered;
            bool x_exterior = myhi(i+1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxlo(i  ,j-1,k  ) == BndryData::covered;
            bool y_exterior = mxlo(i  ,j-1,k  ) == BndryData::not_covered;
            bool z_interior = mxlo(i  ,j  ,k+1) == BndryData::covered;
            bool z_exterior = mxlo(i  ,j  ,k+1) == BndryData::not_covered;
            if (mxlo(i,j,k) != BndryData::covered &&
                (!xlo_domain || !yhi_domain || !zlo_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 3: {
            // xhi & yhi & zlo
            int i = vhi.x+1;
            int j = vhi.y+1;
            int k = vlo.z-1;
            bool x_interior = myhi(i-1,j  ,k  ) == BndryData::covered;
            bool x_exterior = myhi(i-1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j-1,k  ) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j-1,k  ) == BndryData::not_covered;
            bool z_interior = mxhi(i  ,j  ,k+1) == BndryData::covered;
            bool z_exterior = mxhi(i  ,j  ,k+1) == BndryData::not_covered;
            if (mxhi(i,j,k) != BndryData::covered &&
                (!xhi_domain || !yhi_domain || !zlo_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::low, i,j,k, blen.z, vel, mzlo,
                                       bct(Orientation::zlo(), icomp),
                                       bcl(Orientation::zlo(), icomp),
                                       bcvalzlo, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 4: {
            // xlo & ylo & zhi
            int i = vlo.x-1;
            int j = vlo.y-1;
            int k = vhi.z+1;
            bool x_interior = mylo(i+1,j  ,k  ) == BndryData::covered;
            bool x_exterior = mylo(i+1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxlo(i  ,j+1,k  ) == BndryData::covered;
            bool y_exterior = mxlo(i  ,j+1,k  ) == BndryData::not_covered;
            bool z_interior = mxlo(i  ,j  ,k-1) == BndryData::covered;
            bool z_exterior = mxlo(i  ,j  ,k-1) == BndryData::not_covered;
            if (mxlo(i,j,k) != BndryData::covered &&
                (!xlo_domain || !ylo_domain || !zhi_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 5: {
            // xhi & ylo & zhi
            int i = vhi.x+1;
            int j = vlo.y-1;
            int k = vhi.z+1;
            bool x_interior = mylo(i-1,j  ,k  ) == BndryData::covered;
            bool x_exterior = mylo(i-1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j+1,k  ) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j+1,k  ) == BndryData::not_covered;
            bool z_interior = mxhi(i  ,j  ,k-1) == BndryData::covered;
            bool z_exterior = mxhi(i  ,j  ,k-1) == BndryData::not_covered;
            if (mxhi(i,j,k) != BndryData::covered &&
                (!xhi_domain || !ylo_domain || !zhi_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::low, i,j,k, blen.y, vel, mylo,
                                       bct(Orientation::ylo(), icomp),
                                       bcl(Orientation::ylo(), icomp),
                                       bcvalylo, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 6: {
            // xlo & yhi & zhi
            int i = vlo.x-1;
            int j = vhi.y+1;
            int k = vhi.z+1;
            bool x_interior = myhi(i+1,j  ,k  ) == BndryData::covered;
            bool x_exterior = myhi(i+1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxlo(i  ,j-1,k  ) == BndryData::covered;
            bool y_exterior = mxlo(i  ,j-1,k  ) == BndryData::not_covered;
            bool z_interior = mxlo(i  ,j  ,k-1) == BndryData::covered;
            bool z_exterior = mxlo(i  ,j  ,k-1) == BndryData::not_covered;
            if (mxlo(i,j,k) != BndryData::covered &&
                (!xlo_domain || !yhi_domain || !zhi_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::low, i,j,k, blen.x, vel, mxlo,
                                       bct(Orientation::xlo(), icomp),
                                       bcl(Orientation::xlo(), icomp),
                                       bcvalxlo, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        case 7: {
            // xhi & yhi & zhi
            int i = vhi.x+1;
            int j = vhi.y+1;
            int k = vhi.z+1;
            bool x_interior = myhi(i-1,j  ,k  ) == BndryData::covered;
            bool x_exterior = myhi(i-1,j  ,k  ) == BndryData::not_covered;
            bool y_interior = mxhi(i  ,j-1,k  ) == BndryData::covered;
            bool y_exterior = mxhi(i  ,j-1,k  ) == BndryData::not_covered;
            bool z_interior = mxhi(i  ,j  ,k-1) == BndryData::covered;
            bool z_exterior = mxhi(i  ,j  ,k-1) == BndryData::not_covered;
            if (mxhi(i,j,k) != BndryData::covered &&
                (!xhi_domain || !yhi_domain || !zhi_domain)) {
                if ((x_interior && y_interior && z_interior) ||
                    (x_exterior && y_exterior && z_exterior)) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    tmp += vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = Real(1./3.)*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && y_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior && z_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (y_interior && z_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                    Real tmp = vel(i,j,k,icomp);
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                    vel(i,j,k,icomp) = 0.5_rt*(tmp+vel(i,j,k,icomp));
                } else if (x_interior) {
                    mllinop_apply_bc_x(Orientation::high, i,j,k, blen.x, vel, mxhi,
                                       bct(Orientation::xhi(), icomp),
                                       bcl(Orientation::xhi(), icomp),
                                       bcvalxhi, maxorder, dxinv[0], inhomog, icomp);
                } else if (y_interior) {
                    mllinop_apply_bc_y(Orientation::high, i,j,k, blen.y, vel, myhi,
                                       bct(Orientation::yhi(), icomp),
                                       bcl(Orientation::yhi(), icomp),
                                       bcvalyhi, maxorder, dxinv[1], inhomog, icomp);
                } else if (z_interior) {
                    mllinop_apply_bc_z(Orientation::high, i,j,k, blen.z, vel, mzhi,
                                       bct(Orientation::zhi(), icomp),
                                       bcl(Orientation::zhi(), icomp),
                                       bcvalzhi, maxorder, dxinv[2], inhomog, icomp);
                }
            }
            break;
        }
        default: {}
        }
    }
}
#endif

inline
void mltensor_fill_edges (Box const& vbox, // vbox: the valid box
                          Array4<Real> const& vel,
                          Array4<int const> const& mxlo,
                          Array4<int const> const& mylo,
                          Array4<int const> const& mzlo,
                          Array4<int const> const& mxhi,
                          Array4<int const> const& myhi,
                          Array4<int const> const& mzhi,
                          Array4<Real const> const& bcvalxlo,
                          Array4<Real const> const& bcvalylo,
                          Array4<Real const> const& bcvalzlo,
                          Array4<Real const> const& bcvalxhi,
                          Array4<Real const> const& bcvalyhi,
                          Array4<Real const> const& bcvalzhi,
                          Array2D<BoundCond,
                                  0,2*AMREX_SPACEDIM,
                                  0,AMREX_SPACEDIM> const& bct,
                          Array2D<Real,
                                  0,2*AMREX_SPACEDIM,
                                  0,AMREX_SPACEDIM> const& bcl,
                          int inhomog, int maxorder,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Dim3 const& dlo, Dim3 const& dhi) noexcept

{
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    bool xlo_domain = (vlo.x == dlo.x);
    bool ylo_domain = (vlo.y == dlo.y);
    bool zlo_domain = (vlo.z == dlo.z);
    bool xhi_domain = (vhi.x == dhi.x);
    bool yhi_domain = (vhi.y == dhi.y);
    bool zhi_domain = (vhi.z == dhi.z);

    for (int k = vlo.z; k <= vhi.z; ++k) {
        mltensor_fill_edges_xlo_ylo(vlo.x-1, vlo.y-1, k, blen, vel, mxlo, mylo, bcvalxlo, bcvalylo,
                                    bct, bcl, inhomog, maxorder, dxinv, xlo_domain, ylo_domain);
        mltensor_fill_edges_xhi_ylo(vhi.x+1, vlo.y-1, k, blen, vel, mxhi, mylo, bcvalxhi, bcvalylo,
                                    bct, bcl, inhomog, maxorder, dxinv, xhi_domain, ylo_domain);
        mltensor_fill_edges_xlo_yhi(vlo.x-1, vhi.y+1, k, blen, vel, mxlo, myhi, bcvalxlo, bcvalyhi,
                                    bct, bcl, inhomog, maxorder, dxinv, xlo_domain, yhi_domain);
        mltensor_fill_edges_xhi_yhi(vhi.x+1, vhi.y+1, k, blen, vel, mxhi, myhi, bcvalxhi, bcvalyhi,
                                    bct, bcl, inhomog, maxorder, dxinv, xhi_domain, yhi_domain);
    }

    for (int j = vlo.y; j <= vhi.y; ++j) {
        mltensor_fill_edges_xlo_zlo(vlo.x-1, j, vlo.z-1, blen, vel, mxlo, mzlo, bcvalxlo, bcvalzlo,
                                    bct, bcl, inhomog, maxorder, dxinv, xlo_domain, zlo_domain);
        mltensor_fill_edges_xhi_zlo(vhi.x+1, j, vlo.z-1, blen, vel, mxhi, mzlo, bcvalxhi, bcvalzlo,
                                    bct, bcl, inhomog, maxorder, dxinv, xhi_domain, zlo_domain);
        mltensor_fill_edges_xlo_zhi(vlo.x-1, j, vhi.z+1, blen, vel, mxlo, mzhi, bcvalxlo, bcvalzhi,
                                    bct, bcl, inhomog, maxorder, dxinv, xlo_domain, zhi_domain);
        mltensor_fill_edges_xhi_zhi(vhi.x+1, j, vhi.z+1, blen, vel, mxhi, mzhi, bcvalxhi, bcvalzhi,
                                    bct, bcl, inhomog, maxorder, dxinv, xhi_domain, zhi_domain);
    }

    for (int i = vlo.x; i <= vhi.x; ++i) {
        mltensor_fill_edges_ylo_zlo(i, vlo.y-1, vlo.z-1, blen, vel, mylo, mzlo, bcvalylo, bcvalzlo,
                                    bct, bcl, inhomog, maxorder, dxinv, ylo_domain, zlo_domain);
        mltensor_fill_edges_yhi_zlo(i, vhi.y+1, vlo.z-1, blen, vel, myhi, mzlo, bcvalyhi, bcvalzlo,
                                    bct, bcl, inhomog, maxorder, dxinv, yhi_domain, zlo_domain);
        mltensor_fill_edges_ylo_zhi(i, vlo.y-1, vhi.z+1, blen, vel, mylo, mzhi, bcvalylo, bcvalzhi,
                                    bct, bcl, inhomog, maxorder, dxinv, ylo_domain, zhi_domain);
        mltensor_fill_edges_yhi_zhi(i, vhi.y+1, vhi.z+1, blen, vel, myhi, mzhi, bcvalyhi, bcvalzhi,
                                    bct, bcl, inhomog, maxorder, dxinv, yhi_domain, zhi_domain);
    }
}

#ifdef AMREX_USE_GPU
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_edges (int const bid, int const tid, int const bdim,
                          Box const& vbox, // vbox: the valid box
                          Array4<Real> const& vel,
                          Array4<int const> const& mxlo,
                          Array4<int const> const& mylo,
                          Array4<int const> const& mzlo,
                          Array4<int const> const& mxhi,
                          Array4<int const> const& myhi,
                          Array4<int const> const& mzhi,
                          Array4<Real const> const& bcvalxlo,
                          Array4<Real const> const& bcvalylo,
                          Array4<Real const> const& bcvalzlo,
                          Array4<Real const> const& bcvalxhi,
                          Array4<Real const> const& bcvalyhi,
                          Array4<Real const> const& bcvalzhi,
                          Array2D<BoundCond,
                                  0,2*AMREX_SPACEDIM,
                                  0,AMREX_SPACEDIM> const& bct,
                          Array2D<Real,
                                  0,2*AMREX_SPACEDIM,
                                  0,AMREX_SPACEDIM> const& bcl,
                          int inhomog, int maxorder,
                          GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                          Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    bool xlo_domain = (vlo.x == dlo.x);
    bool ylo_domain = (vlo.y == dlo.y);
    bool zlo_domain = (vlo.z == dlo.z);
    bool xhi_domain = (vhi.x == dhi.x);
    bool yhi_domain = (vhi.y == dhi.y);
    bool zhi_domain = (vhi.z == dhi.z);
    if (bid == 0) {
        for (int k = vlo.z + tid; k <= vhi.z; k += bdim) {
            mltensor_fill_edges_xlo_ylo(vlo.x-1, vlo.y-1, k, blen, vel, mxlo, mylo, bcvalxlo, bcvalylo,
                                        bct, bcl, inhomog, maxorder, dxinv, xlo_domain, ylo_domain);
        }
    } else if (bid == 1) {
        for (int k = vlo.z + tid; k <= vhi.z; k += bdim) {
            mltensor_fill_edges_xhi_ylo(vhi.x+1, vlo.y-1, k, blen, vel, mxhi, mylo, bcvalxhi, bcvalylo,
                                        bct, bcl, inhomog, maxorder, dxinv, xhi_domain, ylo_domain);
        }
    } else if (bid == 2) {
        for (int k = vlo.z + tid; k <= vhi.z; k += bdim) {
            mltensor_fill_edges_xlo_yhi(vlo.x-1, vhi.y+1, k, blen, vel, mxlo, myhi, bcvalxlo, bcvalyhi,
                                        bct, bcl, inhomog, maxorder, dxinv, xlo_domain, yhi_domain);
        }
    } else if (bid == 3) {
        for (int k = vlo.z + tid; k <= vhi.z; k += bdim) {
            mltensor_fill_edges_xhi_yhi(vhi.x+1, vhi.y+1, k, blen, vel, mxhi, myhi, bcvalxhi, bcvalyhi,
                                        bct, bcl, inhomog, maxorder, dxinv, xhi_domain, yhi_domain);
        }
    } else if (bid == 4) {
        for (int j = vlo.y + tid; j <= vhi.y; j += bdim) {
            mltensor_fill_edges_xlo_zlo(vlo.x-1, j, vlo.z-1, blen, vel, mxlo, mzlo, bcvalxlo, bcvalzlo,
                                        bct, bcl, inhomog, maxorder, dxinv, xlo_domain, zlo_domain);
        }
    } else if (bid == 5) {
        for (int j = vlo.y + tid; j <= vhi.y; j += bdim) {
            mltensor_fill_edges_xhi_zlo(vhi.x+1, j, vlo.z-1, blen, vel, mxhi, mzlo, bcvalxhi, bcvalzlo,
                                        bct, bcl, inhomog, maxorder, dxinv, xhi_domain, zlo_domain);
        }
    } else if (bid == 6) {
        for (int j = vlo.y + tid; j <= vhi.y; j += bdim) {
            mltensor_fill_edges_xlo_zhi(vlo.x-1, j, vhi.z+1, blen, vel, mxlo, mzhi, bcvalxlo, bcvalzhi,
                                        bct, bcl, inhomog, maxorder, dxinv, xlo_domain, zhi_domain);
        }
    } else if (bid == 7) {
        for (int j = vlo.y + tid; j <= vhi.y; j += bdim) {
            mltensor_fill_edges_xhi_zhi(vhi.x+1, j, vhi.z+1, blen, vel, mxhi, mzhi, bcvalxhi, bcvalzhi,
                                        bct, bcl, inhomog, maxorder, dxinv, xhi_domain, zhi_domain);
        }
    } else if (bid == 8) {
        for (int i = vlo.x + tid; i <= vhi.x; i += bdim) {
            mltensor_fill_edges_ylo_zlo(i, vlo.y-1, vlo.z-1, blen, vel, mylo, mzlo, bcvalylo, bcvalzlo,
                                        bct, bcl, inhomog, maxorder, dxinv, ylo_domain, zlo_domain);
        }
    } else if (bid == 9) {
        for (int i = vlo.x + tid; i <= vhi.x; i += bdim) {
            mltensor_fill_edges_yhi_zlo(i, vhi.y+1, vlo.z-1, blen, vel, myhi, mzlo, bcvalyhi, bcvalzlo,
                                        bct, bcl, inhomog, maxorder, dxinv, yhi_domain, zlo_domain);
        }
    } else if (bid == 10) {
        for (int i = vlo.x + tid; i <= vhi.x; i += bdim) {
            mltensor_fill_edges_ylo_zhi(i, vlo.y-1, vhi.z+1, blen, vel, mylo, mzhi, bcvalylo, bcvalzhi,
                                        bct, bcl, inhomog, maxorder, dxinv, ylo_domain, zhi_domain);
        }
    } else if (bid == 11) {
        for (int i = vlo.x + tid; i <= vhi.x; i += bdim) {
            mltensor_fill_edges_yhi_zhi(i, vhi.y+1, vhi.z+1, blen, vel, myhi, mzhi, bcvalyhi, bcvalzhi,
                                        bct, bcl, inhomog, maxorder, dxinv, yhi_domain, zhi_domain);
        }
    }
}
#endif

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dz_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dzi) noexcept
{
    return (vel(i,j,k+1,n)+vel(i-1,j,k+1,n)-vel(i,j,k-1,n)-vel(i-1,j,k-1,n))*(Real(0.25)*dzi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dz_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dzi) noexcept
{
    return (vel(i,j,k+1,n)+vel(i,j-1,k+1,n)-vel(i,j,k-1,n)-vel(i,j-1,k-1,n))*(Real(0.25)*dzi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_zface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi) noexcept
{
    return (vel(i+1,j,k,n)+vel(i+1,j,k-1,n)-vel(i-1,j,k,n)-vel(i-1,j,k-1,n))*(Real(0.25)*dxi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_zface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi) noexcept
{
    return (vel(i,j+1,k,n)+vel(i,j+1,k-1,n)-vel(i,j-1,k,n)-vel(i,j-1,k-1,n))*(Real(0.25)*dyi);
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi);
                Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi);
                Real dudz = mltensor_dz_on_xface(i,j,k,0,vel,dzi);
                Real dwdz = mltensor_dz_on_xface(i,j,k,2,vel,dzi);
                Real divu = dvdy + dwdz;
                Real xif = kapx(i,j,k);
                Real mun = Real(0.75)*(etax(i,j,k,0)-xif);  // restore the original eta
                Real mut =             etax(i,j,k,1);
                fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,k,1) = -mut*(dudy);
                fx(i,j,k,2) = -mut*(dudz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi);
                Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi);
                Real dvdz = mltensor_dz_on_yface(i,j,k,1,vel,dzi);
                Real dwdz = mltensor_dz_on_yface(i,j,k,2,vel,dzi);
                Real divu = dudx + dwdz;
                Real xif = kapy(i,j,k);
                Real mun = Real(0.75)*(etay(i,j,k,1)-xif);  // restore the original eta
                Real mut =             etay(i,j,k,0);
                fy(i,j,k,0) = -mut*(dvdx);
                fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                fy(i,j,k,2) = -mut*(dvdz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etaz,
                              Array4<Real const> const& kapz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_zface(i,j,k,0,vel,dxi);
                Real dwdx = mltensor_dx_on_zface(i,j,k,2,vel,dxi);
                Real dvdy = mltensor_dy_on_zface(i,j,k,1,vel,dyi);
                Real dwdy = mltensor_dy_on_zface(i,j,k,2,vel,dyi);
                Real divu = dudx + dvdy;
                Real xif = kapz(i,j,k);
                Real mun = Real(0.75)*(etaz(i,j,k,2)-xif);  // restore the original eta
                Real mut =             etaz(i,j,k,0);
                fz(i,j,k,0) = -mut*(dwdx);
                fz(i,j,k,1) = -mut*(dwdy);
                fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dz_on_xface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dzi,
                           Array4<Real const> const& bvxlo, Array4<Real const> const& bvxhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddz;
    if (i == dlo.x) {
        if (bct(Orientation::xlo(),n) == AMREX_LO_DIRICHLET && bvxlo) {
            if (k == dlo.z) {
                ddz = (bvxlo(i-1,j,k  ,n) * Real(-1.5) +
                       bvxlo(i-1,j,k+1,n) * Real(2.) +
                       bvxlo(i-1,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvxlo(i-1,j,k  ,n) * Real(-1.5) +
                        bvxlo(i-1,j,k-1,n) * Real(2.) +
                        bvxlo(i-1,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = (bvxlo(i-1,j,k+1,n)-bvxlo(i-1,j,k-1,n))*(Real(0.5)*dzi);
            }
        } else if (bct(Orientation::xlo(),n) == AMREX_LO_NEUMANN) {
            ddz = (vel(i,j,k+1,n)-vel(i,j,k-1,n))*(Real(0.5)*dzi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else if (i == dhi.x+1) {
        if (bct(Orientation::xhi(),n) == AMREX_LO_DIRICHLET && bvxhi) {
            if (k == dlo.z) {
                ddz = (bvxhi(i,j,k  ,n) * Real(-1.5) +
                       bvxhi(i,j,k+1,n) * Real(2.) +
                       bvxhi(i,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvxhi(i,j,k  ,n) * Real(-1.5) +
                        bvxhi(i,j,k-1,n) * Real(2.) +
                        bvxhi(i,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = (bvxhi(i,j,k+1,n)-bvxhi(i,j,k-1,n))*(Real(0.5)*dzi);
            }
        } else if (bct(Orientation::xhi(),n) == AMREX_LO_NEUMANN) {
            ddz = (vel(i-1,j,k+1,n)-vel(i-1,j,k-1,n))*(Real(0.5)*dzi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else {
        ddz = mltensor_dz_on_xface(i,j,k,n,vel,dzi);
    }
    return ddz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dz_on_yface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dzi,
                           Array4<Real const> const& bvylo, Array4<Real const> const& bvyhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddz;
    if (j == dlo.y) {
        if (bct(Orientation::ylo(),n) == AMREX_LO_DIRICHLET && bvylo) {
            if (k == dlo.z) {
                ddz = (bvylo(i,j-1,k  ,n) * Real(-1.5) +
                       bvylo(i,j-1,k+1,n) * Real(2.) +
                       bvylo(i,j-1,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvylo(i,j-1,k  ,n) * Real(-1.5) +
                        bvylo(i,j-1,k-1,n) * Real(2.) +
                        bvylo(i,j-1,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = (bvylo(i,j-1,k+1,n)-bvylo(i,j-1,k-1,n))*(Real(0.5)*dzi);
            }
        } else if (bct(Orientation::ylo(),n) == AMREX_LO_NEUMANN) {
            ddz = (vel(i,j,k+1,n)-vel(i,j,k-1,n))*(Real(0.5)*dzi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else if (j == dhi.y+1) {
        if (bct(Orientation::yhi(),n) == AMREX_LO_DIRICHLET && bvyhi) {
            if (k == dlo.z) {
                ddz = (bvyhi(i,j,k  ,n) * Real(-1.5) +
                       bvyhi(i,j,k+1,n) * Real(2.) +
                       bvyhi(i,j,k+2,n) * Real(-0.5)) * dzi;
            } else if (k == dhi.z) {
                ddz = -(bvyhi(i,j,k  ,n) * Real(-1.5) +
                        bvyhi(i,j,k-1,n) * Real(2.) +
                        bvyhi(i,j,k-2,n) * Real(-0.5)) * dzi;
            } else {
                ddz = (bvyhi(i,j,k+1,n)-bvyhi(i,j,k-1,n))*(Real(0.5)*dzi);
            }
        } else if (bct(Orientation::yhi(),n) == AMREX_LO_NEUMANN) {
            ddz = (vel(i,j-1,k+1,n)-vel(i,j-1,k-1,n))*(Real(0.5)*dzi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddz = Real(0.);
        }
    } else {
        ddz = mltensor_dz_on_yface(i,j,k,n,vel,dzi);
    }
    return ddz;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dx_on_zface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dxi,
                           Array4<Real const> const& bvzlo, Array4<Real const> const& bvzhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddx;
    if (k == dlo.z) {
        if (bct(Orientation::zlo(),n) == AMREX_LO_DIRICHLET && bvzlo) {
            if (i == dlo.x) {
                ddx = (bvzlo(i  ,j,k-1,n) * Real(-1.5) +
                       bvzlo(i+1,j,k-1,n) * Real(2.) +
                       bvzlo(i+2,j,k-1,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvzlo(i  ,j,k-1,n) * Real(-1.5) +
                        bvzlo(i-1,j,k-1,n) * Real(2.) +
                        bvzlo(i-2,j,k-1,n) * Real(-0.5)) * dxi;
            } else {
                ddx = (bvzlo(i+1,j,k-1,n)-bvzlo(i-1,j,k-1,n))*(Real(0.5)*dxi);
            }
        } else if (bct(Orientation::zlo(),n) == AMREX_LO_NEUMANN) {
            ddx = (vel(i+1,j,k,n)-vel(i-1,j,k,n))*(Real(0.5)*dxi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else if (k == dhi.z+1) {
        if (bct(Orientation::zhi(),n) == AMREX_LO_DIRICHLET && bvzhi) {
            if (i == dlo.x) {
                ddx = (bvzhi(i  ,j,k,n) * Real(-1.5) +
                       bvzhi(i+1,j,k,n) * Real(2.) +
                       bvzhi(i+2,j,k,n) * Real(-0.5)) * dxi;
            } else if (i == dhi.x) {
                ddx = -(bvzhi(i  ,j,k,n) * Real(-1.5) +
                        bvzhi(i-1,j,k,n) * Real(2.) +
                        bvzhi(i-2,j,k,n) * Real(-0.5)) * dxi;
            } else {
                ddx = (bvzhi(i+1,j,k,n)-bvzhi(i-1,j,k,n))*(Real(0.5)*dxi);
            }
        } else if (bct(Orientation::zhi(),n) == AMREX_LO_NEUMANN) {
            ddx = (vel(i+1,j,k-1,n)-vel(i-1,j,k-1,n))*(Real(0.5)*dxi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddx = Real(0.);
        }
    } else {
        ddx = mltensor_dx_on_zface(i,j,k,n,vel,dxi);
    }
    return ddx;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
Real mltensor_dy_on_zface (int i, int j, int k, int n, Array4<Real const> const& vel, Real dyi,
                           Array4<Real const> const& bvzlo, Array4<Real const> const& bvzhi,
                           Array2D<BoundCond,
                                   0,2*AMREX_SPACEDIM,
                                   0,AMREX_SPACEDIM> const& bct,
                           Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    Real ddy;
    if (k == dlo.z) {
        if (bct(Orientation::zlo(),n) == AMREX_LO_DIRICHLET && bvzlo) {
            if (j == dlo.y) {
                ddy = (bvzlo(i,j  ,k-1,n) * Real(-1.5) +
                       bvzlo(i,j+1,k-1,n) * Real(2.) +
                       bvzlo(i,j+2,k-1,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvzlo(i,j  ,k-1,n) * Real(-1.5) +
                        bvzlo(i,j-1,k-1,n) * Real(2.) +
                        bvzlo(i,j-2,k-1,n) * Real(-0.5)) * dyi;
            } else {
                ddy = (bvzlo(i,j+1,k-1,n)-bvzlo(i,j-1,k-1,n))*(Real(0.5)*dyi);
            }
        } else if (bct(Orientation::zlo(),n) == AMREX_LO_NEUMANN) {
            ddy = (vel(i,j+1,k,n)-vel(i,j-1,k,n))*(Real(0.5)*dyi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else if (k == dhi.z+1) {
        if (bct(Orientation::zhi(),n) == AMREX_LO_DIRICHLET && bvzhi) {
            if (j == dlo.y) {
                ddy = (bvzhi(i,j  ,k,n) * Real(-1.5) +
                       bvzhi(i,j+1,k,n) * Real(2.) +
                       bvzhi(i,j+2,k,n) * Real(-0.5)) * dyi;
            } else if (j == dhi.y) {
                ddy = -(bvzhi(i,j  ,k,n) * Real(-1.5) +
                        bvzhi(i,j-1,k,n) * Real(2.) +
                        bvzhi(i,j-2,k,n) * Real(-0.5)) * dyi;
            } else {
                ddy = (bvzhi(i,j+1,k,n)-bvzhi(i,j-1,k,n))*(Real(0.5)*dyi);
            }
        } else if (bct(Orientation::zhi(),n) == AMREX_LO_NEUMANN) {
            ddy = (vel(i,j+1,k-1,n)-vel(i,j-1,k-1,n))*(Real(0.5)*dyi);
        } else { // AMREX_LO_REFLECT_ODD or homogeneous Dirichlet
            ddy = Real(0.);
        }
    } else {
        ddy = mltensor_dy_on_zface(i,j,k,n,vel,dyi);
    }
    return ddy;
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvxlo,
                              Array4<Real const> const& bvxhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dudz = mltensor_dz_on_xface(i,j,k,0,vel,dzi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dwdz = mltensor_dz_on_xface(i,j,k,2,vel,dzi,bvxlo,bvxhi,bct,dlo,dhi);
                Real divu = dvdy + dwdz;
                Real xif = kapx(i,j,k);
                Real mun = Real(0.75)*(etax(i,j,k,0)-xif);  // restore the original eta
                Real mut =             etax(i,j,k,1);
                fx(i,j,k,0) = -mun*(-twoThirds*divu) - xif*divu;
                fx(i,j,k,1) = -mut*(dudy);
                fx(i,j,k,2) = -mut*(dudz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvylo,
                              Array4<Real const> const& bvyhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
                Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
                Real dvdz = mltensor_dz_on_yface(i,j,k,1,vel,dzi,bvylo,bvyhi,bct,dlo,dhi);
                Real dwdz = mltensor_dz_on_yface(i,j,k,2,vel,dzi,bvylo,bvyhi,bct,dlo,dhi);
                Real divu = dudx + dwdz;
                Real xif = kapy(i,j,k);
                Real mun = Real(0.75)*(etay(i,j,k,1)-xif);  // restore the original eta
                Real mut =             etay(i,j,k,0);
                fy(i,j,k,0) = -mut*(dvdx);
                fy(i,j,k,1) = -mun*(-twoThirds*divu) - xif*divu;
                fy(i,j,k,2) = -mut*(dvdz);
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etaz,
                              Array4<Real const> const& kapz,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Array4<Real const> const& bvzlo,
                              Array4<Real const> const& bvzhi,
                              Array2D<BoundCond,
                                      0,2*AMREX_SPACEDIM,
                                      0,AMREX_SPACEDIM> const& bct,
                              Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = Real(2./3.);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_zface(i,j,k,0,vel,dxi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dwdx = mltensor_dx_on_zface(i,j,k,2,vel,dxi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dvdy = mltensor_dy_on_zface(i,j,k,1,vel,dyi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dwdy = mltensor_dy_on_zface(i,j,k,2,vel,dyi,bvzlo,bvzhi,bct,dlo,dhi);
                Real divu = dudx + dvdy;
                Real xif = kapz(i,j,k);
                Real mun = Real(0.75)*(etaz(i,j,k,2)-xif);  // restore the original eta
                Real mut =             etaz(i,j,k,0);
                fz(i,j,k,0) = -mut*(dwdx);
                fz(i,j,k,1) = -mut*(dwdy);
                fz(i,j,k,2) = -mun*(-twoThirds*divu) - xif*divu;
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           Array4<Real const> const& fz,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                           Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const Real dzi = bscalar * dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                Ax(i,j,k,0) += dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                    +          dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                    +          dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0));
                Ax(i,j,k,1) += dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                    +          dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                    +          dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1));
                Ax(i,j,k,2) += dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                    +          dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                    +          dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2));
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_os (Box const& box, Array4<Real> const& Ax,
                              Array4<Real const> const& fx,
                              Array4<Real const> const& fy,
                              Array4<Real const> const& fz,
                              Array4<int const> const& osm,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                              Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const Real dzi = bscalar * dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {
                if (osm(i,j,k) == 0) {
                    Ax(i,j,k,0) = 0.0;
                    Ax(i,j,k,1) = 0.0;
                    Ax(i,j,k,2) = 0.0;
                } else {
                    Ax(i,j,k,0) += dxi*(fx(i+1,j  ,k  ,0) - fx(i,j,k,0))
                        +          dyi*(fy(i  ,j+1,k  ,0) - fy(i,j,k,0))
                        +          dzi*(fz(i  ,j  ,k+1,0) - fz(i,j,k,0));
                    Ax(i,j,k,1) += dxi*(fx(i+1,j  ,k  ,1) - fx(i,j,k,1))
                        +          dyi*(fy(i  ,j+1,k  ,1) - fy(i,j,k,1))
                        +          dzi*(fz(i  ,j  ,k+1,1) - fz(i,j,k,1));
                    Ax(i,j,k,2) += dxi*(fx(i+1,j  ,k  ,2) - fx(i,j,k,2))
                        +          dyi*(fy(i  ,j+1,k  ,2) - fy(i,j,k,2))
                        +          dzi*(fz(i  ,j  ,k+1,2) - fz(i,j,k,2));
                }
            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {

                Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;

                Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi);
                Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi);
                Real dwdy = mltensor_dy_on_xface(i,j,k,2,vel,dyi);

                Real dudz = mltensor_dz_on_xface(i,j,k,0,vel,dzi);
                Real dvdz = mltensor_dz_on_xface(i,j,k,1,vel,dzi);
                Real dwdz = mltensor_dz_on_xface(i,j,k,2,vel,dzi);

                fx(i,j,k,0) = dudx;
                fx(i,j,k,1) = dvdx;
                fx(i,j,k,2) = dwdx;
                fx(i,j,k,3) = dudy;
                fx(i,j,k,4) = dvdy;
                fx(i,j,k,5) = dwdy;
                fx(i,j,k,6) = dudz;
                fx(i,j,k,7) = dvdz;
                fx(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {

                Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi);
                Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi);
                Real dwdx = mltensor_dx_on_yface(i,j,k,2,vel,dxi);

                Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;

                Real dudz = mltensor_dz_on_yface(i,j,k,0,vel,dzi);
                Real dvdz = mltensor_dz_on_yface(i,j,k,1,vel,dzi);
                Real dwdz = mltensor_dz_on_yface(i,j,k,2,vel,dzi);

                fy(i,j,k,0) = dudx;
                fy(i,j,k,1) = dvdx;
                fy(i,j,k,2) = dwdx;
                fy(i,j,k,3) = dudy;
                fy(i,j,k,4) = dvdy;
                fy(i,j,k,5) = dwdy;
                fy(i,j,k,6) = dudz;
                fy(i,j,k,7) = dvdz;
                fy(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                              Array4<Real const> const& vel,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            AMREX_PRAGMA_SIMD
            for (int i = lo.x; i <= hi.x; ++i) {

                Real dudx = mltensor_dx_on_zface(i,j,k,0,vel,dxi);
                Real dvdx = mltensor_dx_on_zface(i,j,k,1,vel,dxi);
                Real dwdx = mltensor_dx_on_zface(i,j,k,2,vel,dxi);

                Real dudy = mltensor_dy_on_zface(i,j,k,0,vel,dyi);
                Real dvdy = mltensor_dy_on_zface(i,j,k,1,vel,dyi);
                Real dwdy = mltensor_dy_on_zface(i,j,k,2,vel,dyi);

                Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;

                fz(i,j,k,0) = dudx;
                fz(i,j,k,1) = dvdx;
                fz(i,j,k,2) = dwdx;
                fz(i,j,k,3) = dudy;
                fz(i,j,k,4) = dvdy;
                fz(i,j,k,5) = dwdy;
                fz(i,j,k,6) = dudz;
                fz(i,j,k,7) = dvdz;
                fz(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fx (Box const& box, Array4<Real> const& fx,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Array4<Real const> const& bvxlo,
                            Array4<Real const> const& bvxhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = (vel(i,j,k,0) - vel(i-1,j,k,0))*dxi;
                Real dvdx = (vel(i,j,k,1) - vel(i-1,j,k,1))*dxi;
                Real dwdx = (vel(i,j,k,2) - vel(i-1,j,k,2))*dxi;
                Real dudy = mltensor_dy_on_xface(i,j,k,0,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dvdy = mltensor_dy_on_xface(i,j,k,1,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dwdy = mltensor_dy_on_xface(i,j,k,2,vel,dyi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dudz = mltensor_dz_on_xface(i,j,k,0,vel,dzi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dvdz = mltensor_dz_on_xface(i,j,k,1,vel,dzi,bvxlo,bvxhi,bct,dlo,dhi);
                Real dwdz = mltensor_dz_on_xface(i,j,k,2,vel,dzi,bvxlo,bvxhi,bct,dlo,dhi);
                fx(i,j,k,0) = dudx;
                fx(i,j,k,1) = dvdx;
                fx(i,j,k,2) = dwdx;
                fx(i,j,k,3) = dudy;
                fx(i,j,k,4) = dvdy;
                fx(i,j,k,5) = dwdy;
                fx(i,j,k,6) = dudz;
                fx(i,j,k,7) = dvdz;
                fx(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fy (Box const& box, Array4<Real> const& fy,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Array4<Real const> const& bvylo,
                            Array4<Real const> const& bvyhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_yface(i,j,k,0,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
                Real dvdx = mltensor_dx_on_yface(i,j,k,1,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
                Real dwdx = mltensor_dx_on_yface(i,j,k,2,vel,dxi,bvylo,bvyhi,bct,dlo,dhi);
                Real dudy = (vel(i,j,k,0) - vel(i,j-1,k,0))*dyi;
                Real dvdy = (vel(i,j,k,1) - vel(i,j-1,k,1))*dyi;
                Real dwdy = (vel(i,j,k,2) - vel(i,j-1,k,2))*dyi;
                Real dudz = mltensor_dz_on_yface(i,j,k,0,vel,dzi,bvylo,bvyhi,bct,dlo,dhi);
                Real dvdz = mltensor_dz_on_yface(i,j,k,1,vel,dzi,bvylo,bvyhi,bct,dlo,dhi);
                Real dwdz = mltensor_dz_on_yface(i,j,k,2,vel,dzi,bvylo,bvyhi,bct,dlo,dhi);
                fy(i,j,k,0) = dudx;
                fy(i,j,k,1) = dvdx;
                fy(i,j,k,2) = dwdx;
                fy(i,j,k,3) = dudy;
                fy(i,j,k,4) = dvdy;
                fy(i,j,k,5) = dwdy;
                fy(i,j,k,6) = dudz;
                fy(i,j,k,7) = dvdz;
                fy(i,j,k,8) = dwdz;

            }
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_vel_grads_fz (Box const& box, Array4<Real> const& fz,
                            Array4<Real const> const& vel,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                            Array4<Real const> const& bvzlo,
                            Array4<Real const> const& bvzhi,
                            Array2D<BoundCond,
                                    0,2*AMREX_SPACEDIM,
                                    0,AMREX_SPACEDIM> const& bct,
                            Dim3 const& dlo, Dim3 const& dhi) noexcept
{
    const Real dxi = dxinv[0];
    const Real dyi = dxinv[1];
    const Real dzi = dxinv[2];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for         (int k = lo.z; k <= hi.z; ++k) {
        for     (int j = lo.y; j <= hi.y; ++j) {
            for (int i = lo.x; i <= hi.x; ++i) {
                Real dudx = mltensor_dx_on_zface(i,j,k,0,vel,dxi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dvdx = mltensor_dx_on_zface(i,j,k,1,vel,dxi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dwdx = mltensor_dx_on_zface(i,j,k,2,vel,dxi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dudy = mltensor_dy_on_zface(i,j,k,0,vel,dyi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dvdy = mltensor_dy_on_zface(i,j,k,1,vel,dyi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dwdy = mltensor_dy_on_zface(i,j,k,2,vel,dyi,bvzlo,bvzhi,bct,dlo,dhi);
                Real dudz = (vel(i,j,k,0) - vel(i,j,k-1,0))*dzi;
                Real dvdz = (vel(i,j,k,1) - vel(i,j,k-1,1))*dzi;
                Real dwdz = (vel(i,j,k,2) - vel(i,j,k-1,2))*dzi;
                fz(i,j,k,0) = dudx;
                fz(i,j,k,1) = dvdx;
                fz(i,j,k,2) = dwdx;
                fz(i,j,k,3) = dudy;
                fz(i,j,k,4) = dvdy;
                fz(i,j,k,5) = dwdy;
                fz(i,j,k,6) = dudz;
                fz(i,j,k,7) = dvdz;
                fz(i,j,k,8) = dwdz;

            }
        }
    }
}

}

#endif
