
#include <AMReX_FAmrCore.H>
#include <AMReX_Print.H>
#include <string>

amrex::FAmrCore::FAmrCore ()
    : amrex::AmrCore()
{
    for (int lev = 0; lev <= maxLevel(); ++lev)
    {
        for (int idim = 0; idim < BL_SPACEDIM; ++idim)
        {
            if (maxGridSize(lev)[idim] % blockingFactor(lev)[idim] != 0 &&
                blockingFactor(lev)[idim] % maxGridSize(lev)[idim] != 0)
            {
                amrex::Abort("On level " + std::to_string(lev) 
                             + " amr.max_grid_size = " + std::to_string(maxGridSize(lev)[idim])
                             + " is not a multiple of amr.blocking_factor = "
                             + std::to_string(blockingFactor(lev)[idim]));
            }
        }
            
        if (blockingFactor(lev) < IntVect{AMREX_D_DECL(8,8,8)})
        {
            if (verbose) {
                amrex::Print() << "amr.blocking_factor < 8 not recommended\n";
            }
        }
    }

    for (int lev = 0; lev < maxLevel(); ++lev)
    {
        const IntVect& rr = refRatio(lev);
        if (!(rr == rr[0]))
        {
            amrex::Abort("On level " + std::to_string(lev) + " amr.ref_ratio is different for different directions."
);
        }
    }

    for (int lev = 0; lev < maxLevel(); ++lev)
    {
        for (int idim = 0; idim < BL_SPACEDIM; ++idim)
        {
            const int rr = refRatio(lev)[idim];
            const int bf = blockingFactor(lev+1)[idim];
            
            if (bf % rr != 0)
            {
                amrex::Abort(" blocking_fact = " + std::to_string(bf)
                             + " is not a multiple of ref_ratio = " + std::to_string(rr));
            }
            
            if (bf / rr < 2)
            {
                amrex::Abort(" blocking_fact = " + std::to_string(bf)
                             + " is too small relative to ref_ratio = " + std::to_string(rr));
            }
        }
    }
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
amrex::FAmrCore::ErrorEst (int lev, TagBoxArray& tags, Real time, int)
{
    if (error_est != nullptr)
    {
	const char   tagval = TagBox::SET;
	const char clearval = TagBox::CLEAR;
	error_est(lev, &tags, time, tagval, clearval);
    }
    else
    {
	amrex::Abort("FAmrCore::error_est is null");
    }
}
