
#include <AMReX_FAmrCore.H>


amrex::FAmrCore::FAmrCore ()
    : amrex::AmrCore()
{
}

amrex::FAmrCore::~FAmrCore ()
{
}

void
amrex::FAmrCore::MakeNewLevelFromScratch (int lev, Real time, const BoxArray& ba,
					  const DistributionMapping& dm)
{
    if (make_new_level_from_scratch != nullptr) {
	make_new_level_from_scratch(lev, time, &ba, &dm);
    } else {
	amrex::Abort("FAmrCore::make_new_level_from_scratch is null");
    }
}

void
amrex::FAmrCore::MakeNewLevelFromCoarse (int lev, Real time, const BoxArray& ba,
					 const DistributionMapping& dm)
{
    if (make_new_level_from_coarse != nullptr) {
	make_new_level_from_coarse(lev, time, &ba, &dm);
    } else {
	amrex::Abort("FAmrCore::make_new_level_from_coarse is null");
    }
}

void
amrex::FAmrCore::RemakeLevel (int lev, Real time, const BoxArray& ba,
			      const DistributionMapping& dm)
{
    if (remake_level != nullptr) {
	remake_level(lev, time, &ba, &dm);
    } else {
	amrex::Abort("FAmrCore::remake_level is null");
    }
}

void
amrex::FAmrCore::ClearLevel (int lev)
{
    if (clear_level != nullptr) {
	clear_level(lev);
    } else {
	amrex::Abort("FAmrCore::clear_level is null");
    }
}

void
amrex::FAmrCore::ErrorEst (int lev, TagBoxArray& tags, Real time, int ngrow)
{
    if (error_est != nullptr) {
	error_est(lev, &tags, time, ngrow);
    } else {
	amrex::Abort("FAmrCore::error_est is null");
    }
}
