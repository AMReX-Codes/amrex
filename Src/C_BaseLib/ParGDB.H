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
	    BoxArray            & ba);

    ParGDB (PArray<Geometry>            & geom, 
	    PArray<DistributionMapping> & dmap,
	    PArray<BoxArray>            & ba,
	    Array<int>                  & rr);

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

    PArray<Geometry>            m_geom;
    PArray<DistributionMapping> m_dmap;
    PArray<BoxArray>            m_ba;
    Array<int>                  m_rr;
    int                         m_nlevels;

};

inline
ParGDB::ParGDB (Geometry            & geom, 
		DistributionMapping & dmap,
		BoxArray            & ba)
    :
    m_nlevels(1)
{
    m_geom.resize(1);
    m_geom.set(0, &geom);

    m_dmap.resize(1);
    m_dmap.set(0, &dmap);

    m_ba.resize(1);
    m_ba.set(0, &ba);
}

inline
ParGDB::ParGDB (PArray<Geometry>            & geom, 
		PArray<DistributionMapping> & dmap,
		PArray<BoxArray>            & ba,
		Array<int>                  & rr)
    :
    m_rr(rr),
    m_nlevels(ba.size())
{
    m_geom.resize(m_nlevels);
    m_dmap.resize(m_nlevels);
    m_ba.resize(m_nlevels);
    for (int lev=0; lev<m_nlevels; ++lev) {
	m_geom.set(lev, &geom[lev]);
	m_dmap.set(lev, &dmap[lev]);
	m_ba.set  (lev, &ba[lev]);
    }
}

inline 
const Geometry& 
ParGDB::Geom (int level) const
{
    return m_geom[level];
}

inline 
const DistributionMapping& 
ParGDB::ParticleDistributionMap (int level) const
{
    return m_dmap[level];
}

inline 
const BoxArray&
ParGDB::ParticleBoxArray (int level) const
{
    return m_ba[level];
}

inline
const BoxArray&
ParGDB::boxArray (int level) const
{
    return m_ba[level];
}

inline 
bool 
ParGDB::LevelDefined (int level) const
{
    return (level < m_nlevels);
}

inline 
int 
ParGDB::finestLevel () const
{
    return m_nlevels-1;
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
    for (int lev = 0; lev < m_nlevels-1; lev++)
       max_ref_ratio = std::max(max_ref_ratio, m_rr[lev]);
    return max_ref_ratio;
}

inline
bool 
ParGDB::ParticlesOnSameGrids (int level) const
{
    return true;
}

#endif


