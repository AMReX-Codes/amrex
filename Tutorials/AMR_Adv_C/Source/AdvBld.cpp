
#include <LevelBld.H>
#include <Adv.H>

class AdvBld
    :
    public LevelBld
{
    virtual void variableSetUp () BL_OVERRIDE;
    virtual void variableCleanUp () BL_OVERRIDE;
    virtual AmrLevel *operator() () BL_OVERRIDE;
    virtual AmrLevel *operator() (Amr&            papa,
                                  int             lev,
                                  const Geometry& level_geom,
                                  const BoxArray& ba,
                                  Real            time) BL_OVERRIDE;
};

AdvBld Adv_bld;

LevelBld*
getLevelBld ()
{
    return &Adv_bld;
}

void
AdvBld::variableSetUp ()
{
    Adv::variableSetUp();
}

void
AdvBld::variableCleanUp ()
{
    Adv::variableCleanUp();
}

AmrLevel*
AdvBld::operator() ()
{
    return new Adv;
}

AmrLevel*
AdvBld::operator() (Amr&            papa,
		    int             lev,
		    const Geometry& level_geom,
		    const BoxArray& ba,
		    Real            time)
{
    return new Adv(papa, lev, level_geom, ba, time);
}
