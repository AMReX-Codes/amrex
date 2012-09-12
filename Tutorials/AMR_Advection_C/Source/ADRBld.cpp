
#include "LevelBld.H"
#include "ADR.H"

class ADRBld
    :
    public LevelBld
{
    virtual void variableSetUp ();
    virtual void variableCleanUp ();
    virtual AmrLevel *operator() ();
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time);
};

ADRBld ADR_bld;

LevelBld*
getLevelBld ()
{
    return &ADR_bld;
}

void
ADRBld::variableSetUp ()
{
    ADR::variableSetUp();
}

void
ADRBld::variableCleanUp ()
{
    ADR::variableCleanUp();
}

AmrLevel*
ADRBld::operator() ()
{
    return new ADR;
}

AmrLevel*
ADRBld::operator() (Amr&            papa,
                       int             lev,
                       const Geometry& level_geom,
                       const BoxArray& ba,
                       Real            time)
{
    return new ADR(papa, lev, level_geom, ba, time);
}
