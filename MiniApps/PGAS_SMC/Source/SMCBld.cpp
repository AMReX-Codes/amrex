
#include "LevelBld.H"
#include "SMC.H"

class SMCBld
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

SMCBld SMC_bld;

LevelBld*
getLevelBld ()
{
    return &SMC_bld;
}

void
SMCBld::variableSetUp ()
{
    SMC::variableSetUp();
}

void
SMCBld::variableCleanUp ()
{
    SMC::variableCleanUp();
}

AmrLevel*
SMCBld::operator() ()
{
    return new SMC;
}

AmrLevel*
SMCBld::operator() (Amr&            papa,
		    int             lev,
		    const Geometry& level_geom,
		    const BoxArray& ba,
		    Real            time)
{
    return new SMC(papa, lev, level_geom, ba, time);
}
