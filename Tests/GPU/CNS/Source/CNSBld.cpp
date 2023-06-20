
#include <AMReX_LevelBld.H>
#include <CNS.H>

using namespace amrex;

class CNSBld
    :
    public LevelBld
{
    virtual void variableSetUp () override;
    virtual void variableCleanUp () override;
    virtual AmrLevel *operator() () override;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  const DistributionMapping& dm,
                                  Real            time) override;
};

CNSBld CNS_bld;

LevelBld*
getLevelBld ()
{
    return &CNS_bld;
}

void
CNSBld::variableSetUp ()
{
    CNS::variableSetUp();
}

void
CNSBld::variableCleanUp ()
{
    CNS::variableCleanUp();
}

AmrLevel*
CNSBld::operator() ()
{
    return new CNS;
}

AmrLevel*
CNSBld::operator() (Amr&            papa,
                    int             lev,
                    const Geometry& level_geom,
                    const BoxArray& ba,
                    const DistributionMapping& dm,
                    Real            time)
{
    return new CNS(papa, lev, level_geom, ba, dm, time);
}
