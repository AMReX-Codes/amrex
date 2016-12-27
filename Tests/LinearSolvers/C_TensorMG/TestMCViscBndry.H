#ifndef _MCVISCBNDRY_H_
#define _MCVISCBNDRY_H_

#include <AMReX_MCInterpBndryData.H>
#include <AMReX_BC_TYPES.H>

class MCViscBndry : public amrex::MCInterpBndryData
{
public:
    MCViscBndry() : amrex::MCInterpBndryData()  {}

    MCViscBndry(const amrex::BoxArray& _grids, 
		const amrex::DistributionMapping& _dmap,
		const amrex::Geometry& geom) :
#if BL_SPACEDIM == 2
        amrex::MCInterpBndryData(_grids,_dmap,4,geom) {};
#elif BL_SPACEDIM == 3
        amrex::MCInterpBndryData(_grids,_dmap,3*(1+3),geom) {};
#endif

    virtual void setBndryConds (const amrex::BCRec& phys_bc,
                                int          ratio,
                                int          comp=0);
    void setHomogValues();
};

#endif
