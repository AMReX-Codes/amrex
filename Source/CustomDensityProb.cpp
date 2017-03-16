#include <PlasmaInjector.H>

#include <iostream>

#include <AMReX.H>

amrex::Real CustomDensityDistribution::getDensity(amrex::Real x,
                                                  amrex::Real y,
                                                  amrex::Real z) const
{
    amrex::Abort("If running with a custom density profile, you must supply a CustomDensityProb.cpp file");
    return 0.0;
}
