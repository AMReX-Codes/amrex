#include <AMReX_ParmParse.H>
#include <AMReX_GpuMemory.H>
#include <AMReX_Vector.H>

#include "AmrLevelAdv.H"

void
AmrLevelAdv::get_tagging_params()
{
    // Use ParmParse to get number of levels from input file
    amrex::ParmParse pp("tagging");
    pp.query("max_phierr_lev", max_phierr_lev);
    pp.query("max_phigrad_lev", max_phigrad_lev);

    // Set default values for the error thresholds, then read from input file
    if (max_phierr_lev  != -1) {
        phierr.resize(max_phierr_lev, 1.0e+20);
        pp.queryarr("phierr",  phierr);
    }
    if (max_phigrad_lev != -1) {
        phigrad.resize(max_phigrad_lev, 1.0e+20);
        pp.queryarr("phigrad", phigrad);
    }
}
