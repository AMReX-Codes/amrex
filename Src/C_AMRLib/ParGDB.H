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

    ParGDB (const Geometry            & geom, 
	    const DistributionMapping & dmap,
	    const BoxArray            & ba)
    {
	m_geom.resize(1);
	m_dmap.resize(1);
	m_ba.resize(1);
	m_geom.set(0, &geom);  
	m_dmap.set(0, &dmap);  
	m_ba.set(0, &ba);  	
    }

    ParGDB (const PArray<Geometry>            & geom,
	    const PArray<DistributionMapping> & dmap,
	    const PArray<BoxArray>            & ba)
    {
	BL_ASSERT(geom.size() == dmap.size() && geom.size() == ba.size());
	int N = geom.size();
	m_geom.resize(N);
	m_dmap.resize(N);
	m_ba.resize(N);
	for (int i = 0; i < N; ++i) {
	    m_geom.set(i, &geom[i]);  
	    m_dmap.set(i, &dmap[i]);  
	    m_ba.set(i, &ba[i]);  
	}
    }

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
    // What do we do with ref ratio?
    // Should we pass it in the multi-level constructor?
    virtual IntVect refRatio (int level) const;
    virtual int MaxRefRatio (int level) const;

    virtual bool ParticlesOnSameGrids (int level) const;

protected:

    PArray<Geometry>            m_geom;
    PArray<DistributionMapping> m_dmap;
    PArray<BoxArray>            m_ba;

};

inline 
const Geometry& 
ParGDB::Geom (int level) const
{
    BL_ASSERT(level < m_geom.size());
    return m_geom[0];
}

inline 
const DistributionMapping& 
ParGDB::ParticleDistributionMap (int level) const
{
    BL_ASSERT(level < m_dmap.size());
    return m_dmap[0];
}

inline 
const BoxArray&
ParGDB::ParticleBoxArray (int level) const
{
    BL_ASSERT(level < m_ba.size());
    return m_ba[0];
}

inline
const BoxArray&
ParGDB::boxArray (int level) const
{
    BL_ASSERT(level < m_ba.size());
    return m_ba[0];
}

inline 
bool 
ParGDB::LevelDefined (int level) const
{
    return level < m_geom.size();
}

inline 
int 
ParGDB::finestLevel () const
{
    return m_geom.size()-1;
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
