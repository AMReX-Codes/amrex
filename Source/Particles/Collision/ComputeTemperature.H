/* Copyright 2019-2020 Andrew Myers, Yinjian Zhao
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARTICLES_COLLISION_COMPUTE_TEMPERATURE_H_
#define WARPX_PARTICLES_COLLISION_COMPUTE_TEMPERATURE_H_

#include "Utils/WarpXConst.H"


template <typename T_index, typename T_R>
AMREX_GPU_HOST_DEVICE
T_R ComputeTemperature (
    T_index const Is, T_index const Ie, T_index const *I,
    T_R const *ux, T_R const *uy, T_R const *uz, T_R const m )
{

    T_R constexpr inv_c2 = T_R(1.0) / ( PhysConst::c * PhysConst::c );

    int N = Ie - Is;
    if ( N == 0 ) { return T_R(0.0); }

    T_R vx = T_R(0.0);    T_R vy = T_R(0.0);
    T_R vz = T_R(0.0);    T_R vs = T_R(0.0);
    T_R gm = T_R(0.0);    T_R us = T_R(0.0);

    for (int i = Is; i < Ie; ++i)
    {
        us = ( ux[ I[i] ] * ux[ I[i] ] +
               uy[ I[i] ] * uy[ I[i] ] +
               uz[ I[i] ] * uz[ I[i] ] );
        gm = std::sqrt( T_R(1.0) + us*inv_c2 );
        vx += ux[ I[i] ] / gm;
        vy += uy[ I[i] ] / gm;
        vz += uz[ I[i] ] / gm;
        vs += us / gm / gm;
    }

    vx = vx / N;    vy = vy / N;
    vz = vz / N;    vs = vs / N;

    return m/T_R(3.0)*(vs-(vx*vx+vy*vy+vz*vz));
}

#endif // WARPX_PARTICLES_COLLISION_COMPUTE_TEMPERATURE_H_
