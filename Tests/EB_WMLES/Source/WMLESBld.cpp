
#include <AMReX_LevelBld.H>
#include <WMLES.H>

using namespace amrex;

class WMLESBld
    :
    public LevelBld
{
    void variableSetUp () override;
    void variableCleanUp () override;
    AmrLevel *operator() () override;
    AmrLevel *operator() (Amr&            papa,
                          int             lev,
                          const Geometry& level_geom,
                          const BoxArray& ba,
                          const DistributionMapping& dm,
                          Real            time) override;
};

WMLESBld wmles_bld;

LevelBld*
getLevelBld ()
{
    return &wmles_bld;
}

void
WMLESBld::variableSetUp ()
{
    WMLES::variableSetUp();
}

void
WMLESBld::variableCleanUp ()
{
    WMLES::variableCleanUp();
}

AmrLevel*
WMLESBld::operator() ()
{
    return new WMLES;
}

AmrLevel*
WMLESBld::operator() (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    const DistributionMapping& dm,
                    Real            time)
{
    return new WMLES(papa, lev, level_geom, ba, dm, time);
}
