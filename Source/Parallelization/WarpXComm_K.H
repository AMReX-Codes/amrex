/* Copyright 2019 Axel Huebl, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_COMM_K_H_
#define WARPX_COMM_K_H_

#include <AMReX_FArrayBox.H>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Bxa,
                            amrex::Array4<amrex::Real const> const& Bxf,
                            amrex::Array4<amrex::Real const> const& Bxc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    Real const wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real const owx = 1.0_rt - wx;
    Bxa(j,k,l) = owx * Bxc(jg,kg,lg) + wx * Bxc(jg+1,kg,lg) + Bxf(j,k,l);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Bya,
                            amrex::Array4<amrex::Real const> const& Byf,
                            amrex::Array4<amrex::Real const> const& Byc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real const wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real const owy = 1.0_rt - wy;
    Bya(j,k,l) = owy * Byc(jg,kg,lg) + wy * Byc(jg,kg+1,lg) + Byf(j,k,l);
#else
    Bya(j,k,l) = Byc(jg,kg,lg) + Byf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_bfield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Bza,
                            amrex::Array4<amrex::Real const> const& Bzf,
                            amrex::Array4<amrex::Real const> const& Bzc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    // Note that for 2d, l=0, because the amrex convention is used here.

#if (AMREX_SPACEDIM == 3)
    Real const wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real const owz = 1.0_rt - wz;
    Bza(j,k,l) = owz * Bzc(jg,kg,lg) + owz * Bzc(jg,kg,lg+1) + Bzf(j,k,l);
#else
    Real const wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real const owy = 1.0_rt - wy;
    Bza(j,k,l) = owy * Bzc(jg,kg,lg) + owy * Bzc(jg,kg+1,lg) + Bzf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_x (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Exa,
                            amrex::Array4<amrex::Real const> const& Exf,
                            amrex::Array4<amrex::Real const> const& Exc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    Real const wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real const owy = 1.0_rt - wy;

#if (AMREX_SPACEDIM == 3)
    Real const wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real const owz = 1.0_rt - wz;
    Exa(j,k,l) = owy * owz * Exc(jg  ,kg  ,lg  )
        +         wy * owz * Exc(jg  ,kg+1,lg  )
        +        owy *  wz * Exc(jg  ,kg  ,lg+1)
        +         wy *  wz * Exc(jg  ,kg+1,lg+1)
        + Exf(j,k,l);
#else
    Exa(j,k,l) = owy * Exc(jg,kg,lg) + wy * Exc(jg,kg+1,lg) + Exf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_y (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Eya,
                            amrex::Array4<amrex::Real const> const& Eyf,
                            amrex::Array4<amrex::Real const> const& Eyc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    Real const wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real const owx = 1.0_rt - wx;

#if (AMREX_SPACEDIM == 3)
    Real const wz = (l == lg*2) ? 0.0_rt : 0.5_rt;
    Real const owz = 1.0_rt - wz;
    Eya(j,k,l) = owx * owz * Eyc(jg  ,kg  ,lg  )
        +         wx * owz * Eyc(jg+1,kg  ,lg  )
        +        owx *  wz * Eyc(jg  ,kg  ,lg+1)
        +         wx *  wz * Eyc(jg+1,kg  ,lg+1)
        + Eyf(j,k,l);
#else
    Real const wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real const owy = 1.0_rt - wy;
    Eya(j,k,l) = owx * owy * Eyc(jg  ,kg  ,lg)
        +         wx * owy * Eyc(jg+1,kg  ,lg)
        +        owx *  wy * Eyc(jg  ,kg+1,lg)
        +         wx *  wy * Eyc(jg+1,kg+1,lg)
        + Eyf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_efield_z (int j, int k, int l,
                            amrex::Array4<amrex::Real      > const& Eza,
                            amrex::Array4<amrex::Real const> const& Ezf,
                            amrex::Array4<amrex::Real const> const& Ezc)
{
    using namespace amrex;

    int const lg = amrex::coarsen(l,2);
    int const kg = amrex::coarsen(k,2);
    int const jg = amrex::coarsen(j,2);

    Real const wx = (j == jg*2) ? 0.0_rt : 0.5_rt;
    Real const owx = 1.0_rt - wx;

#if (AMREX_SPACEDIM == 3)
    Real wy = (k == kg*2) ? 0.0_rt : 0.5_rt;
    Real owy = 1.0_rt - wy;
    Eza(j,k,l) = owx * owy * Ezc(jg  ,kg  ,lg  )
        +         wx * owy * Ezc(jg+1,kg  ,lg  )
        +        owx *  wy * Ezc(jg  ,kg+1,lg  )
        +         wx *  wy * Ezc(jg+1,kg+1,lg  )
        + Ezf(j,k,l);
#else
    Eza(j,k,l) = owx * Ezc(jg,kg,lg) + wx * Ezc(jg+1,kg,lg) + Ezf(j,k,l);
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf,
                               amrex::Array4<amrex::Real const> const& Bxc,
                               amrex::Array4<amrex::Real const> const& Bxg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bxg(jg  ,kg  ,0)
        +     owx *  wy * Bxg(jg  ,kg+1,0)
        +      wx * owy * Bxg(jg+1,kg  ,0)
        +      wx *  wy * Bxg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    Real bc = owx * owy * Bxc(jg  ,kg  ,0)
        +     owx *  wy * Bxc(jg  ,kg-1,0)
        +      wx * owy * Bxc(jg+1,kg  ,0)
        +      wx *  wy * Bxc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5*(Bxf(j,k-1,0) + Bxf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Bxg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bxg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bxg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Bxg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Bxg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bxg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Bxg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Bxg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    wz = 0.5-wz;  owz = 1.0-wz;
    Real bc = owx * owy * owz * Bxc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bxc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bxc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Bxc(jg+1,kg-1,lg  )
        +     owx * owy *  wz * Bxc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Bxc(jg+1,kg  ,lg-1)
        +     owx *  wy *  wz * Bxc(jg  ,kg-1,lg-1)
        +      wx *  wy *  wz * Bxc(jg+1,kg-1,lg-1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25*(Bxf(j,k-1,l-1) + Bxf(j,k,l-1) + Bxf(j,k-1,l) + Bxf(j,k,l));
#endif

    Bxa(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf,
                               amrex::Array4<amrex::Real const> const& Byc,
                               amrex::Array4<amrex::Real const> const& Byg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Byg(jg  ,kg  ,0)
        +     owx *  wy * Byg(jg  ,kg+1,0)
        +      wx * owy * Byg(jg+1,kg  ,0)
        +      wx *  wy * Byg(jg+1,kg+1,0);

    // interp from coarse stagged (cell-centered for By) to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    wy = 0.5-wy;  owy = 1.0-wy;
    Real bc = owx * owy * Byc(jg  ,kg  ,0)
        +     owx *  wy * Byc(jg  ,kg-1,0)
        +      wx * owy * Byc(jg-1,kg  ,0)
        +      wx *  wy * Byc(jg-1,kg-1,0);

    // interp form fine stagger (cell-centered for By) to fine nodal
    Real bf = 0.25*(Byf(j,k,0) + Byf(j-1,k,0) + Byf(j,k-1,0) + Byf(j-1,k-1,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Byg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Byg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Byg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Byg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Byg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Byg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Byg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Byg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    wz = 0.5-wz;  owz = 1.0-wz;
    Real bc = owx * owy * owz * Byc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Byc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Byc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Byc(jg-1,kg+1,lg  )
        +     owx * owy *  wz * Byc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Byc(jg-1,kg  ,lg-1)
        +     owx *  wy *  wz * Byc(jg  ,kg+1,lg-1)
        +      wx *  wy *  wz * Byc(jg-1,kg+1,lg-1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25*(Byf(j-1,k,l-1) + Byf(j,k,l-1) + Byf(j-1,k,l) + Byf(j,k,l));

#endif

    Bya(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf,
                               amrex::Array4<amrex::Real const> const& Bzc,
                               amrex::Array4<amrex::Real const> const& Bzg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * Bzg(jg  ,kg  ,0)
        +     owx *  wy * Bzg(jg  ,kg+1,0)
        +      wx * owy * Bzg(jg+1,kg  ,0)
        +      wx *  wy * Bzg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    Real bc = owx * owy * Bzc(jg  ,kg  ,0)
        +     owx *  wy * Bzc(jg  ,kg+1,0)
        +      wx * owy * Bzc(jg-1,kg  ,0)
        +      wx *  wy * Bzc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real bf = 0.5*(Bzf(j-1,k,0) + Bzf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real bg = owx * owy * owz * Bzg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bzg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Bzg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Bzg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Bzg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bzg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Bzg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Bzg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    wy = 0.5-wy;  owy = 1.0-wy;
    Real bc = owx * owy * owz * Bzc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Bzc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Bzc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Bzc(jg-1,kg-1,lg  )
        +     owx * owy *  wz * Bzc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Bzc(jg-1,kg  ,lg+1)
        +     owx *  wy *  wz * Bzc(jg  ,kg-1,lg+1)
        +      wx *  wy *  wz * Bzc(jg-1,kg-1,lg+1);

    // interp from fine stagged to fine nodal
    Real bf = 0.25*(Bzf(j-1,k-1,l) + Bzf(j,k-1,l) + Bzf(j-1,k,l) + Bzf(j,k,l));

#endif

    Bza(j,k,l) = bg + (bf-bc);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bxa,
                               amrex::Array4<amrex::Real const> const& Bxf)
{
#if (AMREX_SPACEDIM == 2)
    Bxa(j,k,0) = 0.5*(Bxf(j,k-1,0) + Bxf(j,k,0));
#else
    Bxa(j,k,l) = 0.25*(Bxf(j,k-1,l-1) + Bxf(j,k,l-1) + Bxf(j,k-1,l) + Bxf(j,k,l));
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bya,
                               amrex::Array4<amrex::Real const> const& Byf)
{
#if (AMREX_SPACEDIM == 2)
    Bya(j,k,0) = 0.25*(Byf(j,k,0) + Byf(j-1,k,0) + Byf(j,k-1,0) + Byf(j-1,k-1,0));
#else
    Bya(j,k,l) = 0.25*(Byf(j-1,k,l-1) + Byf(j,k,l-1) + Byf(j-1,k,l) + Byf(j,k,l));
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_bfield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Bza,
                               amrex::Array4<amrex::Real const> const& Bzf)
{
#if (AMREX_SPACEDIM == 2)
    Bza(j,k,0) = 0.5*(Bzf(j-1,k,0) + Bzf(j,k,0));
#else
    Bza(j,k,l) = 0.25*(Bzf(j-1,k-1,l) + Bzf(j,k-1,l) + Bzf(j-1,k,l) + Bzf(j,k,l));
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf,
                               amrex::Array4<amrex::Real const> const& Exc,
                               amrex::Array4<amrex::Real const> const& Exg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Exg(jg  ,kg  ,0)
        +     owx *  wy * Exg(jg  ,kg+1,0)
        +      wx * owy * Exg(jg+1,kg  ,0)
        +      wx *  wy * Exg(jg+1,kg+1,0);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    Real ec = owx * owy * Exc(jg  ,kg  ,0)
        +     owx *  wy * Exc(jg  ,kg+1,0)
        +      wx * owy * Exc(jg-1,kg  ,0)
        +      wx *  wy * Exc(jg-1,kg+1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Exf(j-1,k,0) + Exf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Exg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Exg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Exg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Exg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Exg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Exg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Exg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Exg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wx = 0.5-wx;  owx = 1.0-wx;
    Real ec = owx * owy * owz * Exc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Exc(jg-1,kg  ,lg  )
        +     owx *  wy * owz * Exc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Exc(jg-1,kg+1,lg  )
        +     owx * owy *  wz * Exc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Exc(jg-1,kg  ,lg+1)
        +     owx *  wy *  wz * Exc(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Exc(jg-1,kg+1,lg+1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Exf(j-1,k,l) + Exf(j,k,l));

#endif

    Exa(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf,
                               amrex::Array4<amrex::Real const> const& Eyc,
                               amrex::Array4<amrex::Real const> const& Eyg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal and coarse staggered to fine nodal
    Real eg = owx * owy * (Eyg(jg  ,kg  ,0) + Eyc(jg  ,kg  ,0))
        +     owx *  wy * (Eyg(jg  ,kg+1,0) + Eyc(jg  ,kg+1,0))
        +      wx * owy * (Eyg(jg+1,kg  ,0) + Eyc(jg+1,kg  ,0))
        +      wx *  wy * (Eyg(jg+1,kg+1,0) + Eyc(jg+1,kg+1,0));
    Real ec = 0.0;

    // interp from fine staggered to fine nodal
    Real ef = Eyf(j,k,0);

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Eyg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Eyg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Eyg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Eyg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Eyg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Eyg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Eyg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Eyg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    Real ec = owx * owy * owz * Eyc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Eyc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Eyc(jg  ,kg-1,lg  )
        +      wx *  wy * owz * Eyc(jg+1,kg-1,lg  )
        +     owx * owy *  wz * Eyc(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Eyc(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Eyc(jg  ,kg-1,lg+1)
        +      wx *  wy *  wz * Eyc(jg+1,kg-1,lg+1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Eyf(j,k-1,l) + Eyf(j,k,l));

#endif

    Eya(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf,
                               amrex::Array4<amrex::Real const> const& Ezc,
                               amrex::Array4<amrex::Real const> const& Ezg)
{
    using namespace amrex;

    int jg = amrex::coarsen(j,2);
    Real wx = (j == jg*2) ? 0.0 : 0.5;
    Real owx = 1.0-wx;

    int kg = amrex::coarsen(k,2);
    Real wy = (k == kg*2) ? 0.0 : 0.5;
    Real owy = 1.0-wy;

#if (AMREX_SPACEDIM == 2)

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * Ezg(jg  ,kg  ,0)
        +     owx *  wy * Ezg(jg  ,kg+1,0)
        +      wx * owy * Ezg(jg+1,kg  ,0)
        +      wx *  wy * Ezg(jg+1,kg+1,0);

    // interp from coarse stagged to fine nodal
    wy = 0.5-wy;  owy = 1.0-wy;
    Real ec = owx * owy * Ezc(jg  ,kg  ,0)
        +     owx *  wy * Ezc(jg  ,kg-1,0)
        +      wx * owy * Ezc(jg+1,kg  ,0)
        +      wx *  wy * Ezc(jg+1,kg-1,0);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Ezf(j,k-1,0) + Ezf(j,k,0));

#else

    int lg = amrex::coarsen(l,2);
    Real wz = (l == lg*2) ? 0.0 : 0.5;
    Real owz = 1.0-wz;

    // interp from coarse nodal to fine nodal
    Real eg = owx * owy * owz * Ezg(jg  ,kg  ,lg  )
        +      wx * owy * owz * Ezg(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Ezg(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Ezg(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Ezg(jg  ,kg  ,lg+1)
        +      wx * owy *  wz * Ezg(jg+1,kg  ,lg+1)
        +     owx *  wy *  wz * Ezg(jg  ,kg+1,lg+1)
        +      wx *  wy *  wz * Ezg(jg+1,kg+1,lg+1);

    // interp from coarse staggered to fine nodal
    wz = 0.5-wz;  owz = 1.0-wz;
    Real ec = owx * owy * owz * Ezc(jg  ,kg  ,lg  )
        +      wx * owy * owz * Ezc(jg+1,kg  ,lg  )
        +     owx *  wy * owz * Ezc(jg  ,kg+1,lg  )
        +      wx *  wy * owz * Ezc(jg+1,kg+1,lg  )
        +     owx * owy *  wz * Ezc(jg  ,kg  ,lg-1)
        +      wx * owy *  wz * Ezc(jg+1,kg  ,lg-1)
        +     owx *  wy *  wz * Ezc(jg  ,kg+1,lg-1)
        +      wx *  wy *  wz * Ezc(jg+1,kg+1,lg-1);

    // interp from fine staggered to fine nodal
    Real ef = 0.5*(Ezf(j,k,l-1) + Ezf(j,k,l));

#endif

    Eza(j,k,l) = eg + (ef-ec);
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_x (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Exa,
                               amrex::Array4<amrex::Real const> const& Exf)
{
    Exa(j,k,l) = 0.5*(Exf(j-1,k,l) + Exf(j,k,l));
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_y (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eya,
                               amrex::Array4<amrex::Real const> const& Eyf)
{
#if (AMREX_SPACEDIM == 2)
    Eya(j,k,0) = Eyf(j,k,0);
#else
    Eya(j,k,l) = 0.5*(Eyf(j,k-1,l) + Eyf(j,k,l));
#endif
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void warpx_interp_nd_efield_z (int j, int k, int l,
                               amrex::Array4<amrex::Real> const& Eza,
                               amrex::Array4<amrex::Real const> const& Ezf)
{
#if (AMREX_SPACEDIM == 2)
    Eza(j,k,0) = 0.5*(Ezf(j,k-1,0) + Ezf(j,k,0));
#else
    Eza(j,k,l) = 0.5*(Ezf(j,k,l-1) + Ezf(j,k,l));
#endif
}

#endif
