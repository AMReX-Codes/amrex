#ifndef _AmrCoreAdv_F_H_
#define _AmrCoreAdv_F_H_

#include <AMReX_BLFort.H>
#include <AMReX_Box.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>

using namespace amrex;

    // these are C++ function definitions converted from Fortran 
/*
    AMREX_GPU_DEVICE
    void initdata(Box const& bx, FArrayBox& phifab, 
                  GeometryData const & geom);
*/
// ==== GET_FACE_VELOCITY =====
/*
    AMREX_GPU_DEVICE
    void get_face_velocity_psi(Box const& bx,
                               const Real time,
                               FArrayBox& psifab,
                               GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_x(Box const& bx,
                             FArrayBox& xvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_y(Box const& bx,
                             FArrayBox& yvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);

    AMREX_GPU_DEVICE
    void get_face_velocity_z(Box const& bx,
                             FArrayBox& zvelfab,
                             const FArrayBox& psifab,
                             GeometryData const& geom);
*/
// ==== ADVECT FUNCTIONS =====
/*
    AMREX_GPU_DEVICE
    void conservative(Box const& bx,
                      const FArrayBox& statein,
                      FArrayBox& stateout,
                      AMREX_D_DECL(FArrayBox& fx,
                                   FArrayBox& fy,
                                   FArrayBox& fz),
                      const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_scale_x(Box const& bx,
                      FArrayBox& fx,
                      const Real& dt,
                      const GpuArray<Real, AMREX_SPACEDIM>& dx);

    AMREX_GPU_DEVICE
    void flux_scale_y(Box const& bx,
                      FArrayBox& fy,
                      const Real& dt,
                      const GpuArray<Real, AMREX_SPACEDIM>& dx);

    AMREX_GPU_DEVICE
    void flux_scale_z(Box const& bx,
                      FArrayBox& fz,
                      const Real& dt,
                      const GpuArray<Real, AMREX_SPACEDIM>& dx);
*/

// ==== COMPUTE FLUX FUNCTIONS =====
/*
    AMREX_GPU_DEVICE
    void flux_x(Box const& bx,
                const FArrayBox& state,
                const FArrayBox& velx,
                FArrayBox& phix,
                const FArrayBox& slope,
                const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_y(Box const& bx,
                const FArrayBox& state,
                const FArrayBox& vely,
                FArrayBox& phiy,
                const FArrayBox& slope,
                const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_z(Box const& bx,
                const FArrayBox& state,
                const FArrayBox& velz,
                FArrayBox& phiz,
                const FArrayBox& slope,
                const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_xy(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phix_y,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_xz(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phix_z,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_yx(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phiy_x,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_yz(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phi_yz,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_zx(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phi_zx,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void flux_zy(Box const& bx,
                 AMREX_D_DECL(const FArrayBox& velx,
                              const FArrayBox& vely,
                              const FArrayBox& velz),
                 AMREX_D_DECL(const FArrayBox& phix,
                              const FArrayBox& phiy,
                              const FArrayBox& phiz),
                 FArrayBox& phi_zy,
                 const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void combine_flux_x(Box const& bx,
                        const FArrayBox& velx,
                        const FArrayBox& vely,
                        const FArrayBox& velz,
                        FArrayBox& phix,
                        const FArrayBox& phiy_z,
                        const FArrayBox& phiz_y,
                        FArrayBox& flxx,
                        const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void combine_flux_y(Box const& bx,
                        const FArrayBox& velx,
                        const FArrayBox& vely,
                        const FArrayBox& velz,
                        FArrayBox& phiy,
                        const FArrayBox& phix_z,
                        const FArrayBox& phiz_x,
                        FArrayBox& flxy,
                        const GpuArray<Real, AMREX_SPACEDIM>& dtdx);

    AMREX_GPU_DEVICE
    void combine_flux_z(Box const& bx,
                        const FArrayBox& velx,
                        const FArrayBox& vely,
                        const FArrayBox& velz,
                        FArrayBox& phiz,
                        const FArrayBox& phix_y,
                        const FArrayBox& phiy_x,
                        FArrayBox& flxz,
                        const GpuArray<Real, AMREX_SPACEDIM>& dtdx);
*/
/*
    void state_error(TagBox& tag,
		     const FArrayBox& state,
		     const int tagval, const int clearval,
                     Box const& bx,
		     const amrex_real* phierr);
*/

// ==== SLOPE FUNCTIONS =====

    AMREX_GPU_DEVICE
    void slopex2(Box const& bx,
                const FArrayBox &qfab,
                FArrayBox &dqfab);

    AMREX_GPU_DEVICE
    void slopex4(Box const& bx,
                 const FArrayBox &qfab,
                 const FArrayBox &dqfab,
                 FArrayBox &dq4fab);

    AMREX_GPU_DEVICE
    void slopey2(Box const& bx,
                const FArrayBox &qfab,
                FArrayBox &dqfab);

    AMREX_GPU_DEVICE
    void slopey4(Box const& bx,
                 const FArrayBox &qfab,
                 const FArrayBox &dqfab,
                 FArrayBox &dq4fab);

    AMREX_GPU_DEVICE
    void slopez2(Box const& bx,
                const FArrayBox &qfab,
                FArrayBox &dqfab);

    AMREX_GPU_DEVICE
    void slopez4(Box const& bx,
                 const FArrayBox &qfab,
                 const FArrayBox &dqfab,
                 FArrayBox &dq4fab);


#ifdef __cplusplus
extern "C" 
{
#endif

    void state_error(int* tag, const int* tag_lo, const int* tag_hi,
		     const BL_FORT_FAB_ARG_3D(state),
		     const int* tagval, const int* clearval,
		     const int* lo, const int* hi,
		     const amrex_real* dx, const amrex_real* problo,
		     const amrex_real* time, const amrex_real* phierr);


    void phifill (BL_FORT_FAB_ARG_3D(state),
		  const int* dlo, const int* dhi,
		  const amrex_real* dx, const amrex_real* glo, 
		  const amrex_real* time, const int* bc);

#ifdef __cplusplus
}
#endif

#endif
