#ifndef AMREX_ALGOIM_H_
#define AMREX_ALGOIM_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>

/**
 * \file AMReX_algoim.H
 *
 * Wrapper routines for the Algoim EB integration library.
 */

namespace amrex::algoim {

//! Compute volume integrals of monomials over cut cells (using \p nghost ghost cells, but
//! effectively clamped to the component-wise minimum of \p nghost and intg.nGrowVect()).
void compute_integrals (MultiFab& intg, int nghost = 100);
//! Overload that accepts an IntVect grow vector.
void compute_integrals (MultiFab& intg, IntVect nghost);

//! Compute surface integrals over EB faces (using \p nghost ghost cells, but
//! effectively clamped to the component-wise minimum of \p nghost and intg.nGrowVect()).
void compute_surface_integrals (MultiFab& sintg, int nghost = 100);
//! Overload accepting anisotropic grow vector.
void compute_surface_integrals (MultiFab& sintg, IntVect nghost);

/// \cond DOXYGEN_IGNORE

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

static constexpr int i_B_x        =  0;
static constexpr int i_B_y        =  1;
static constexpr int i_B_z        =  2;
static constexpr int i_B_x_y      =  3;
static constexpr int i_B_x_z      =  4;
static constexpr int i_B_y_z      =  5;
static constexpr int i_B_xyz      =  6;
static constexpr int numSurfIntgs =  7;

/// \endcond

}

#endif
