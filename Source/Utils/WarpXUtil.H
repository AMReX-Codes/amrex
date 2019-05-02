#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_MultiFab.H>

void ReadBoostedFrameParameters(amrex::Real& gamma_boost, amrex::Real& beta_boost,
                                amrex::Vector<int>& boost_direction);

void ConvertLabParamsToBoost();

void NullifyMF(amrex::MultiFab& mf, int lev, amrex::Real zmin, 
               amrex::Real zmax);
