#ifndef _ParGDB_H
#define _ParGDB_H

#include <Array.H>
#include <Geometry.H>
#include <MultiFab.H>
#include <BL_CXX11.H>
				      
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

    virtual void SetParticleBoxArray (int level,
				      const DistributionMapping& new_dm,
				      const BoxArray           & new_ba) = 0;

    virtual bool LevelDefined (int level) const = 0;
    virtual int finestLevel () const = 0;
    virtual int initialBaLevels () const = 0;

    virtual int subCycle () const = 0;
    virtual IntVect refRatio (int level) const = 0;
    virtual int MaxRefRatio (int level) const = 0;

    bool OnSameGrids (int level, const MultiFab& mf) const;
};

//
// we use this for non-Amr particle code
//
class ParGDB
    : public ParGDBBase
{
public:

    ParGDB () { ; }

    ParGDB (const Geometry            & geom, 
	    const DistributionMapping & dmap,
	    const BoxArray            & ba);

    ParGDB (const Array<Geometry>            & geom, 
	    const Array<DistributionMapping> & dmap,
	    const Array<BoxArray>            & ba,
	    const Array<int>                 & rr);

    virtual ~ParGDB () {;}

    virtual const Geometry& Geom (int level) const BL_OVERRIDE;
    virtual const DistributionMapping& ParticleDistributionMap 
                                             (int level) const BL_OVERRIDE;
    virtual const BoxArray& ParticleBoxArray (int level) const BL_OVERRIDE;
    virtual const BoxArray&         boxArray (int level) const BL_OVERRIDE;

    virtual void SetParticleBoxArray (int level,
				      const DistributionMapping& new_dm,
				      const BoxArray           & new_ba) BL_OVERRIDE;

    virtual bool LevelDefined (int level) const BL_OVERRIDE;
    virtual int finestLevel () const BL_OVERRIDE;
    virtual int initialBaLevels () const BL_OVERRIDE;

    virtual int subCycle () const BL_OVERRIDE;
    virtual IntVect refRatio (int level) const BL_OVERRIDE;
    virtual int MaxRefRatio (int level) const BL_OVERRIDE;

protected:

    Array<Geometry>            m_geom;
    Array<DistributionMapping> m_dmap;
    Array<BoxArray>            m_ba;
    Array<int>                 m_rr;
    int                        m_nlevels;

};

inline
bool
ParGDBBase::OnSameGrids (int level, const MultiFab& mf) const
{
    return (mf.DistributionMap() == ParticleDistributionMap(level) &&
	    mf.boxArray().CellEqual(ParticleBoxArray(level)));
}

inline
ParGDB::ParGDB (const Geometry            & geom, 
		const DistributionMapping & dmap,
		const BoxArray            & ba)
    :
    m_geom(1,geom),
    m_dmap(1,dmap),
    m_ba(1,ba),
    m_nlevels(1)
{ }

inline
ParGDB::ParGDB (const Array<Geometry>            & geom, 
		const Array<DistributionMapping> & dmap,
		const Array<BoxArray>            & ba,
		const Array<int>                 & rr)
    :
    m_geom(geom),
    m_dmap(dmap),
    m_ba(ba),
    m_rr(rr),
    m_nlevels(ba.size())
{ }

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
void 
ParGDB::SetParticleBoxArray (int level,
			     const DistributionMapping& new_dm,
			     const BoxArray           & new_ba)
{
    BL_ASSERT(level < m_nlevels);
    m_dmap[level] = new_dm;
    m_ba[level] = new_ba;
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

#endif


