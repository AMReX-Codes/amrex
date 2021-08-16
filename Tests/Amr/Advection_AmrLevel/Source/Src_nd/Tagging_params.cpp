#include <AMReX_ParmParse.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Vector.H>

#include "AmrLevelAdv.H"

int          AmrLevelAdv::max_phierr_lev  = -1;
int          AmrLevelAdv::max_phigrad_lev = -1;
amrex::Real* AmrLevelAdv::phierr;
amrex::Real* AmrLevelAdv::phigrad;

void
AmrLevelAdv::get_tagging_params()
{
    // Use ParmParse to get number of levels from input file
    int maxlev_in;
    amrex::ParmParse pp("tagging");
    pp.query("max_phierr_lev", max_phierr_lev);
    pp.query("max_phigrad_lev", max_phigrad_lev);
    if (max_phigrad_lev == -1) { max_phigrad_lev = max_phierr_lev; }

    // Allocate arrays to hold error thresholds
    phierr  = (amrex::Real*)malloc(max_phierr_lev*sizeof(amrex::Real));
    phigrad = (amrex::Real*)malloc(max_phigrad_lev*sizeof(amrex::Real));

    // Set default values for the error thresholds
    for (int i = 0; i < max_phierr_lev; ++i)  { phierr[i]  = 1.0e+20; }
    for (int i = 0; i < max_phigrad_lev; ++i) { phigrad[i] = 1.0e+20; }

    // Read error thresholds from input file
    amrex::Vector<amrex::Real> phierr_in(max_phierr_lev);
    amrex::Vector<amrex::Real> phigrad_in(max_phigrad_lev);
    pp.queryarr("phierr",  phierr_in);
    pp.queryarr("phigrad", phigrad_in);
    for (int i = 0; i < max_phierr_lev; ++i)   { phierr[i]  = phierr_in[i]; }
    for (int i = 0; i < max_phigrad_lev; ++i)  { phigrad[i] = phigrad_in[i]; }
}
