#include <AMReX_REAL.H>
#include <AMReX_Vector.H>

void ReadBoostedFrameParameters(amrex::Real& gamma_boost, amrex::Real& beta_boost,
                                amrex::Vector<int>& boost_direction);

void ConvertLabParamsToBoost();
