#include <AMReX_FArrayBox.H>
#include <AMReX_Array.H>

#include "AmrLevelAdv.H"
#include "Kernels.H"

AMREX_GPU_HOST
void
AmrLevelAdv::advect (const amrex::Real /*time*/,
                     const amrex::Box& bx,
                     amrex::GpuArray<amrex::Box,BL_SPACEDIM> nbx,
                     const amrex::FArrayBox& statein,
                     amrex::FArrayBox& stateout,
                     AMREX_D_DECL(const amrex::FArrayBox& vx,
                                  const amrex::FArrayBox& vy,
                                  const amrex::FArrayBox& vz),
                     AMREX_D_DECL(amrex::FArrayBox& fx,
                                  amrex::FArrayBox& fy,
                                  amrex::FArrayBox& fz),
                     amrex::GpuArray<amrex::Real,BL_SPACEDIM> dx,
                     const amrex::Real dt)
{
    using namespace amrex;

    const Box& gbx = amrex::grow(bx, 1);

    // Set up dz
#if (AMREX_SPACEDIM == 3)
    Real dz = dx[2];
#else
    Real dz = 1.0;
#endif

    // Set up dt/dx ratios
    Real dtdx  = dt/dx[0];
    Real dtdy  = dt/dx[1];
#if (AMREX_SPACEDIM == 3)
    Real dtdz  = dt/dx[2];
#endif
    Real hdtdx = Real(0.5)*dtdx;
    Real hdtdy = Real(0.5)*dtdy;
#if (AMREX_SPACEDIM == 3)
    Real hdtdz = Real(0.5)*dtdz;
    Real tdtdx = Real(1.0/3.0)*dtdx;
    Real tdtdy = Real(1.0/3.0)*dtdy;
    Real tdtdz = Real(1.0/3.0)*dtdz;
#endif

    // Set up temporary fab with elixir
    FArrayBox tmpfab;
    int ntmpcomps = (AMREX_SPACEDIM == 2) ? 6 : 14;
    tmpfab.resize(amrex::grow(bx,2),ntmpcomps);
    Elixir tmpeli = tmpfab.elixir();
    int itmp = 0;

    // State in and out
    Array4<Real const> phi_in = statein.array();
    Array4<Real>      phi_out = stateout.array();

    // Face velocities
    Array4<Real const> velx = vx.array();
    Array4<Real const> vely = vy.array();
#if (AMREX_SPACEDIM == 3)
    Array4<Real const> velz = vz.array();
#endif

    // Fluxes
    AMREX_D_TERM(Array4<Real> tfluxx = tmpfab.array(itmp++);,
                 Array4<Real> tfluxy = tmpfab.array(itmp++);,
                 Array4<Real> tfluxz = tmpfab.array(itmp++));

    // Slopes
    Array4<Real> slope2 = tmpfab.array(itmp++);
    Array4<Real const> slope2_c = slope2;
    Array4<Real> slope4 = tmpfab.array(itmp++);
    Array4<Real const> slope4_c = slope4;

    // Compute longitudinal fluxes
    // ===========================

    // x -------------------------
    Array4<Real> phix = tmpfab.array(itmp++);
    Array4<Real const> phix_c = phix;

    // Fromm slope
    amrex::launch(amrex::grow(gbx,Direction::x,1),
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopex2(tbx, slope2, phi_in);
    });

    // Limited fourth order slope
    amrex::launch(gbx,
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopex4(tbx, slope4, phi_in, slope2_c);
    });

    // Longitudinal flux
    Box b = gbx;
    amrex::ParallelFor(b.grow(Direction::x,-1).surroundingNodes(Direction::x),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_x(i, j, k, phix, phi_in, velx, slope4_c, hdtdx);
    });

    // y -------------------------
    Array4<Real> phiy = tmpfab.array(itmp++);
    Array4<Real const> phiy_c = phiy;

    // Fromm slope
    amrex::launch(amrex::grow(gbx,Direction::y,1),
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopey2(tbx, slope2, phi_in);
    });

    // Limited fourth order slope
    amrex::launch(gbx,
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopey4(tbx, slope4, phi_in, slope2_c);
    });

    // Longitudinal flux
    b = gbx;
    amrex::ParallelFor(b.grow(Direction::y,-1).surroundingNodes(Direction::y),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_y(i, j, k, phiy, phi_in, vely, slope4_c, hdtdy);
    });

#if (AMREX_SPACEDIM == 3)
    // z -------------------------
    Array4<Real> phiz = tmpfab.array(itmp++);
    Array4<Real const> phiz_c = phiz;

    // Fromm slope
    amrex::launch(amrex::grow(gbx,Direction::z,1),
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopez2(tbx, slope2, phi_in);
    });

    // Limited fourth order slope
    amrex::launch(gbx,
    [=] AMREX_GPU_DEVICE (const Box& tbx)
    {
        slopez4(tbx, slope4, phi_in, slope2_c);
    });

    // Longitudinal flux
    b = gbx;
    amrex::ParallelFor(b.grow(Direction::z,-1).surroundingNodes(Direction::z),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_z(i, j, k, phiz, phi_in, velz, slope4_c, hdtdz);
    });

    // compute transverse fluxes (3D only)
    // ===================================

    // xy --------------------
    Array4<Real> phix_y = tmpfab.array(itmp++);
    Array4<Real const> phix_y_c = phix_y;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::z,1).surroundingNodes(Direction::x),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_xy(i, j, k, phix_y,
                velx, vely,
                phix_c, phiy_c,
                tdtdy);
    });

    // xz --------------------
    Array4<Real> phix_z = tmpfab.array(itmp++);
    Array4<Real const> phix_z_c = phix_z;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::y,1).surroundingNodes(Direction::x),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_xz(i, j, k, phix_z,
                velx, velz,
                phix_c, phiz_c,
                tdtdz);
    });

    // yx --------------------
    Array4<Real> phiy_x = tmpfab.array(itmp++);
    Array4<Real const> phiy_x_c = phiy_x;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::z,1).surroundingNodes(Direction::y),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_yx(i, j, k, phiy_x,
                velx, vely,
                phix_c, phiy_c,
                tdtdx);
    });

    // yz --------------------
    Array4<Real> phiy_z = tmpfab.array(itmp++);
    Array4<Real const> phiy_z_c = phiy_z;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::x,1).surroundingNodes(Direction::y),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_yz(i, j, k, phiy_z,
                vely, velz,
                phiy_c, phiz_c,
                tdtdz);
    });

    // zx --------------------
    Array4<Real> phiz_x = tmpfab.array(itmp++);
    Array4<Real const> phiz_x_c = phiz_x;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::y,1).surroundingNodes(Direction::z),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_zx(i, j, k, phiz_x,
                velx, velz,
                phix_c, phiz_c,
                tdtdx);
    });

    // zy --------------------
    Array4<Real> phiz_y = tmpfab.array(itmp++);
    Array4<Real const> phiz_y_c = phiz_y;

    b = bx;
    amrex::ParallelFor(b.grow(Direction::x,1).surroundingNodes(Direction::z),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        flux_zy(i, j, k, phiz_y,
                vely, velz,
                phiy_c, phiz_c,
                tdtdy);
    });
#endif

    // final edge states
    // ===========================
    amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::x),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        create_flux_x(i, j, k, tfluxx,
                      AMREX_D_DECL(velx,vely,velz),
#if (AMREX_SPACEDIM == 3)
                      phix_c, phiy_z_c, phiz_y_c,
                      hdtdy, hdtdz);
#else
                      phix_c, phiy_c,
                      hdtdy);
#endif
    });

    amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::y),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        create_flux_y(i, j, k, tfluxy,
                      AMREX_D_DECL(velx,vely,velz),
#if (AMREX_SPACEDIM == 3)
                      phiy_c, phix_z_c, phiz_x_c,
                      hdtdx, hdtdz);
#else
                      phiy_c, phix_c,
                      hdtdx);
#endif
    });

#if (AMREX_SPACEDIM == 3)
    amrex::ParallelFor(amrex::surroundingNodes(bx,Direction::z),
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        create_flux_z(i, j, k, tfluxz,
                      velx, vely, velz,
                      phiz_c, phix_y_c, phiy_x_c,
                      hdtdx, hdtdy);
    });
#endif

    AMREX_ASSERT(itmp == ntmpcomps);

    // compute new state (stateout) and scale fluxes based on face area.
    // ===========================

    AMREX_D_TERM(Array4<Real const> tfluxx_c = tfluxx;,
                 Array4<Real const> tfluxy_c = tfluxy;,
                 Array4<Real const> tfluxz_c = tfluxz);

    // Do a conservative update
    amrex::ParallelFor(bx,
    [=] AMREX_GPU_DEVICE (int i, int j, int k)
    {
        conservative(i, j, k,
                     phi_out, phi_in,
                     AMREX_D_DECL(tfluxx_c,tfluxy_c,tfluxz_c),
                     AMREX_D_DECL(dtdx,dtdy,dtdz));
    });

    // Scale by face area in order to correctly reflux
    amrex::ParallelFor(
        AMREX_D_DECL(amrex::surroundingNodes(bx,Direction::x),
                     amrex::surroundingNodes(bx,Direction::y),
                     amrex::surroundingNodes(bx,Direction::z)),
        AMREX_D_DECL([=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         tfluxx(i,j,k) *= dt*dx[1]*dz;
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         tfluxy(i,j,k) *= dt*dx[0]*dz;
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         tfluxz(i,j,k) *= dt*dx[0]*dx[1];
                     }));

    // Copy fluxes into Flux MultiFab
    AMREX_D_TERM(Array4<Real> fluxx = fx.array();,
                 Array4<Real> fluxy = fy.array();,
                 Array4<Real> fluxz = fz.array());
    amrex::ParallelFor(
        AMREX_D_DECL(nbx[0],
                     nbx[1],
                     nbx[2]),
        AMREX_D_DECL([=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         fluxx(i,j,k) = tfluxx_c(i,j,k);
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         fluxy(i,j,k) = tfluxy_c(i,j,k);
                     },
                     [=] AMREX_GPU_DEVICE (int i, int j, int k)
                     {
                         fluxz(i,j,k) = tfluxz_c(i,j,k);
                     }));
}
