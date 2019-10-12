#ifndef AMREX_ALGOIM_H_
#define AMREX_ALGOIM_H_

#include <AMReX_MultiFab.H>

namespace amrex { namespace algoim {

void compute_integrals (MultiFab& intg, int nghost = 100);
void compute_integrals (MultiFab& intg, IntVect nghost);

static constexpr int i_S_x     =  0;
static constexpr int i_S_y     =  1;
static constexpr int i_S_z     =  2;
static constexpr int i_S_x2    =  3;
static constexpr int i_S_y2    =  4;
static constexpr int i_S_z2    =  5;
static constexpr int i_S_x_y   =  6;
static constexpr int i_S_x_z   =  7;
static constexpr int i_S_y_z   =  8;
static constexpr int i_S_x2_y  =  9;
static constexpr int i_S_x2_z  = 10;
static constexpr int i_S_x_y2  = 11;
static constexpr int i_S_y2_z  = 12;
static constexpr int i_S_x_z2  = 13;
static constexpr int i_S_y_z2  = 14;
static constexpr int i_S_x2_y2 = 15;
static constexpr int i_S_x2_z2 = 16;
static constexpr int i_S_y2_z2 = 17;
static constexpr int i_S_xyz   = 18;
static constexpr int numIntgs  = 19;

}}

#endif
