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

    ParGDB (PArray<Geometry>            & geom, 
	     Array<DistributionMapping> & dmap,
	    PArray<BoxArray>            & ba,
	    int                           rr[])
	: m_geom(&geom), m_dmap(&dmap), m_ba(&ba)
	{ m_rr[0] = rr[0];
	  m_finestLevel = m_geom->size()-1;}

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

    PArray<Geometry>            *m_geom;
     Array<DistributionMapping> *m_dmap;
    PArray<BoxArray>            *m_ba;
    int                          m_rr[];
    int                          m_finestLevel;

};

inline 
const Geometry& 
ParGDB::Geom (int level) const
{
    return (*m_geom)[level];
}

inline 
const DistributionMapping& 
ParGDB::ParticleDistributionMap (int level) const
{
    return (*m_dmap)[level];
}

inline 
const BoxArray&
ParGDB::ParticleBoxArray (int level) const
{
    return (*m_ba)[level];
}

inline
const BoxArray&
ParGDB::boxArray (int level) const
{
    return (*m_ba)[level];
}

inline 
bool 
ParGDB::LevelDefined (int level) const
{
    return (level <= m_finestLevel);
}

inline 
int 
ParGDB::finestLevel () const
{
    return m_finestLevel;
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
    return m_rr[level]*IntVect::TheUnitVector();
}

inline 
int 
ParGDB::MaxRefRatio (int level) const
{
    int max_ref_ratio = 0;
    for (int i = 0; i < m_finestLevel; i++)
       max_ref_ratio = std::max(max_ref_ratio, m_rr[level]);
    return max_ref_ratio;
}

inline
bool 
ParGDB::ParticlesOnSameGrids (int level) const
{
    return true;
}

#endif


