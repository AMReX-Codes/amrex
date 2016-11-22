
#include <LevelBld.H>
#include <Adv.H>

class AdvBld
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
                                  Real            time) override;
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
