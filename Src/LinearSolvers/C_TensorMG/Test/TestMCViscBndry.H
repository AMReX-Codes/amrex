#ifndef _MCVISCBNDRY_H_
#define _MCVISCBNDRY_H_

#include <MCInterpBndryData.H>
#include <BC_TYPES.H>

class MCViscBndry : public MCInterpBndryData
{
public:
    MCViscBndry() : MCInterpBndryData()  {}

    MCViscBndry(const BoxArray& _grids, const Geometry& geom) :
#if BL_SPACEDIM == 2
        MCInterpBndryData(_grids,4,geom) {};
#elif BL_SPACEDIM == 3
        MCInterpBndryData(_grids,3*(1+3),geom) {};
#endif

    virtual void setBndryConds (const BCRec& phys_bc,
                                int          ratio,
                                int          comp=0);
    void setHomogValues();
};

#endif
