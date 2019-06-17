#ifndef AMREX_ParGDB_H_
#define AMREX_ParGDB_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_MultiFab.H>

namespace amrex {

class ParGDBBase
{
public:

    ParGDBBase () {;}
    virtual ~ParGDBBase () {;}

    virtual const Geometry& Geom (int level) const = 0;
    virtual const DistributionMapping& ParticleDistributionMap
                                             (int level) const = 0;
    virtual const DistributionMapping&         DistributionMap
                                             (int level) const = 0;
    virtual const BoxArray& ParticleBoxArray (int level) const = 0;
    virtual const BoxArray&         boxArray (int level) const = 0;

    virtual void SetParticleBoxArray (int level, const BoxArray& new_ba) = 0;
    virtual void SetParticleDistributionMap (int level, const DistributionMapping& new_dm) = 0;

    virtual void ClearParticleBoxArray (int level) = 0;
    virtual void ClearParticleDistributionMap (int level) = 0;

    virtual bool LevelDefined (int level) const = 0;
    virtual int finestLevel () const = 0;
    virtual int maxLevel () const = 0;

    virtual IntVect refRatio (int level) const = 0;
    virtual int MaxRefRatio (int level) const = 0;

    bool OnSameGrids (int level, const MultiFab& mf) const;
};

//
//! we use this for non-Amr particle code
class ParGDB
    : public ParGDBBase
{
public:

    ParGDB () { ; }

    ParGDB (const Geometry            & geom,
	    const DistributionMapping & dmap,
	    const BoxArray            & ba);

    ParGDB (const Vector<Geometry>            & geom,
	    const Vector<DistributionMapping> & dmap,
	    const Vector<BoxArray>            & ba,
	    const Vector<int>                 & rr);

    virtual ~ParGDB () {;}

    virtual const Geometry& Geom (int level) const override;
    virtual const DistributionMapping& ParticleDistributionMap
                                             (int level) const override;
    virtual const DistributionMapping&         DistributionMap
                                             (int level) const override;
    virtual const BoxArray& ParticleBoxArray (int level) const override;
    virtual const BoxArray&         boxArray (int level) const override;

    virtual void SetParticleBoxArray (int level, const BoxArray& new_ba) override;
    virtual void SetParticleDistributionMap (int level,	const DistributionMapping& new_dm) override;

    virtual void ClearParticleBoxArray (int level) override;
    virtual void ClearParticleDistributionMap (int level) override;

    virtual bool LevelDefined (int level) const override;
    virtual int finestLevel () const override;
    virtual int maxLevel () const override;

    virtual IntVect refRatio (int level) const override;
    virtual int MaxRefRatio (int level) const override;

protected:

    Vector<Geometry>            m_geom;
    Vector<DistributionMapping> m_dmap;
    Vector<BoxArray>            m_ba;
    Vector<int>                 m_rr;
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
ParGDB::ParGDB (const Vector<Geometry>            & geom,
		const Vector<DistributionMapping> & dmap,
		const Vector<BoxArray>            & ba,
		const Vector<int>                 & rr)
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
const DistributionMapping&
ParGDB::DistributionMap (int level) const
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
ParGDB::SetParticleBoxArray (int level, const BoxArray& new_ba)
{
    BL_ASSERT(level < m_nlevels);
    m_ba[level] = new_ba;
}

inline
void
ParGDB::SetParticleDistributionMap (int level, const DistributionMapping& new_dm)
{
    BL_ASSERT(level < m_nlevels);
    m_dmap[level] = new_dm;
}

inline
void
ParGDB::ClearParticleBoxArray (int level)
{
    BL_ASSERT(level < m_nlevels);
    m_ba[level] = BoxArray();
}

inline
void
ParGDB::ClearParticleDistributionMap (int level)
{
    BL_ASSERT(level < m_nlevels);
    m_dmap[level] = DistributionMapping();
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
ParGDB::maxLevel () const
{
    return m_nlevels-1;
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

}

#endif
