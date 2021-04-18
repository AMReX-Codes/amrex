
#include <AMReX_LevelBld.H>
#include <AmrLevelAdv.H>

using namespace amrex;

class LevelBldAdv
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

LevelBldAdv Adv_bld;

LevelBld*
getLevelBld ()
{
    return &Adv_bld;
}

void
LevelBldAdv::variableSetUp ()
{
    AmrLevelAdv::variableSetUp();
}

void
LevelBldAdv::variableCleanUp ()
{
    AmrLevelAdv::variableCleanUp();
}

AmrLevel*
LevelBldAdv::operator() ()
{
    return new AmrLevelAdv;
}

AmrLevel*
LevelBldAdv::operator() (Amr&            papa,
                         int             lev,
                         const Geometry& level_geom,
                         const BoxArray& ba,
                         const DistributionMapping& dm,
                         Real            time)
{
    return new AmrLevelAdv(papa, lev, level_geom, ba, dm, time);
}
