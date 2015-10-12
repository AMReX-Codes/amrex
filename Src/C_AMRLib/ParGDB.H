#ifndef _ParGDB_H
#define _ParGDB_H

#include <PArray.H>
#include <Geometry.H>
#include <MultiFab.H>

class ParGDBBase
{
public:

    ParGDBBase () {;}
    virtual ~ParGDBBase () {;}

    virtual const Geometry& Geom (int level) const = 0;
    virtual const DistributionMapping& ParticleDistributionMap 
                                             (int level) const = 0;
    virtual const BoxArray& ParticleBoxArray (int level) const = 0;
    virtual const BoxArray&         boxArray (int level) const = 0;

    virtual bool LevelDefined (int level) const = 0;
    virtual int finestLevel () const = 0;
    virtual int initialBaLevels () const = 0;

    virtual int subCycle () const = 0;
    virtual IntVect refRatio (int level) const = 0;
    virtual int MaxRefRatio (int level) const = 0;

    virtual bool ParticlesOnSameGrids (int level) const = 0;
};

//
// we use this for non-Amr particle code
//
class ParGDB
    : public ParGDBBase
{
public:

    ParGDB () { ; }

    ParGDB (Geometry            & geom, 
	    DistributionMapping & dmap,
	    BoxArray            & ba)
	: m_geom(&geom), m_dmap(&dmap), m_ba(&ba)
	{ ; }

    virtual ~ParGDB () {;}

    virtual const Geometry& Geom (int level) const;
    virtual const DistributionMapping& ParticleDistributionMap 
                                             (int level) const;
    virtual const BoxArray& ParticleBoxArray (int level) const;
    virtual const BoxArray&         boxArray (int level) const;

    virtual bool LevelDefined (int level) const;
    virtual int finestLevel () const;
    virtual int initialBaLevels () const;

    virtual int subCycle () const;
    virtual IntVect refRatio (int level) const;
    virtual int MaxRefRatio (int level) const;

    virtual bool ParticlesOnSameGrids (int level) const;

protected:

    Geometry            *m_geom;
    DistributionMapping *m_dmap;
    BoxArray            *m_ba;

};

inline 
const Geometry& 
ParGDB::Geom (int level) const
{
    BL_ASSERT(level == 0);
    return *m_geom;
}

inline 
const DistributionMapping& 
ParGDB::ParticleDistributionMap (int level) const
{
    BL_ASSERT(level == 0);
    return *m_dmap;
}

inline 
const BoxArray&
ParGDB::ParticleBoxArray (int level) const
{
    BL_ASSERT(level == 0);
    return *m_ba;
}

inline
const BoxArray&
ParGDB::boxArray (int level) const
{
    BL_ASSERT(level == 0);
    return *m_ba;
}

inline 
bool 
ParGDB::LevelDefined (int level) const
{
    return level == 0;
}

inline 
int 
ParGDB::finestLevel () const
{
    return 0;
}
 
inline 
int 
ParGDB::initialBaLevels () const
{
    return 0;
}

inline 
int 
ParGDB::subCycle () const
{
    return false;
}

inline 
IntVect 
ParGDB::refRatio (int level) const
{
    return IntVect::TheZeroVector();
}

inline 
int 
ParGDB::MaxRefRatio (int level) const
{
    return 0;
}

inline
bool 
ParGDB::ParticlesOnSameGrids (int level) const
{
    return true;
}

#endif


