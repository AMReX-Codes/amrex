#ifndef AMREX_AmrParGDB_H
#define AMREX_AmrParGDB_H

#include <AMReX_ParGDB.H>
#include <AMReX_AmrCore.H>

namespace amrex {

class AmrParGDB
    : public ParGDBBase
{
public:

    AmrParGDB (AmrCore* amr) noexcept
	: m_amrcore(amr),
	  m_dmap(amr->maxLevel()+1),
	  m_ba(amr->maxLevel()+1)
	{ }

    virtual ~AmrParGDB () {;}
    
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
    AmrCore* m_amrcore;
    Vector<DistributionMapping> m_dmap;
    Vector<BoxArray>            m_ba;
};

inline 
const Geometry& 
AmrParGDB::Geom (int level) const
{
    return m_amrcore->Geom(level);
}

inline 
const DistributionMapping& 
AmrParGDB::ParticleDistributionMap (int level) const
{
    if (m_dmap[level].empty()) {
	return m_amrcore->DistributionMap(level);
    } else {
	return m_dmap[level];
    }
}

inline 
const DistributionMapping& 
AmrParGDB::DistributionMap (int level) const
{
    return m_amrcore->DistributionMap(level);
}

inline 
const BoxArray&
AmrParGDB::ParticleBoxArray (int level) const
{
    if (m_ba[level].empty()) {
	return m_amrcore->boxArray(level);
    } else {
	return m_ba[level];
    }
}

inline
const BoxArray&
AmrParGDB::boxArray (int level) const
{
    return m_amrcore->boxArray(level);
}

inline
void AmrParGDB::SetParticleBoxArray (int level, const BoxArray& new_ba)
{
    m_ba[level] = new_ba;
}

inline
void AmrParGDB::SetParticleDistributionMap (int level, const DistributionMapping& new_dmap)
{
    m_dmap[level] = new_dmap;
}

inline
void AmrParGDB::ClearParticleBoxArray (int level)
{
    m_ba[level] = BoxArray();
}

inline
void AmrParGDB::ClearParticleDistributionMap (int level)
{
    m_dmap[level] = DistributionMapping();
}

inline 
bool 
AmrParGDB::LevelDefined (int level) const
{
    return m_amrcore->LevelDefined(level);
}

inline 
int 
AmrParGDB::finestLevel () const
{
    return m_amrcore->finestLevel();
}

inline 
int 
AmrParGDB::maxLevel () const
{
    return m_amrcore->maxLevel();
}
 
inline 
IntVect 
AmrParGDB::refRatio (int level) const
{
    return m_amrcore->refRatio(level);
}

inline 
int 
AmrParGDB::MaxRefRatio (int level) const
{
    return m_amrcore->MaxRefRatio(level);
}

}

#endif
