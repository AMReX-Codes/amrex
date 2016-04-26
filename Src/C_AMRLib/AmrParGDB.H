#ifndef _AmrParGDB_H
#define _AmrParGDB_H

#include <ParGDB.H>
#include <Amr.H>

class AmrParGDB
    : public ParGDBBase
{
public:

    AmrParGDB (Amr* amr)
	: m_amr(amr) { BL_ASSERT(amr != 0); }

    virtual ~AmrParGDB () {;}
    
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
    Amr* m_amr;
};

inline 
const Geometry& 
AmrParGDB::Geom (int level) const
{
    return m_amr->Geom(level);
}

inline 
const DistributionMapping& 
AmrParGDB::ParticleDistributionMap (int level) const
{
    return m_amr->getLevel(level).ParticleDistributionMap();
}

inline 
const BoxArray&
AmrParGDB::ParticleBoxArray (int level) const
{
    return m_amr->ParticleBoxArray(level);
}

inline
const BoxArray&
AmrParGDB::boxArray (int level) const
{
    return m_amr->boxArray(level);
}

inline
void AmrParGDB::SetParticleBoxArray (int level,
				     const DistributionMapping& new_dm,
				     const BoxArray           & new_ba)
{
    BoxLib::Abort("AmrParGDB::SetParticleBoxArray: not implemented yet");
}

inline 
bool 
AmrParGDB::LevelDefined (int level) const
{
    return m_amr->getAmrLevels().defined(level);
}

inline 
int 
AmrParGDB::finestLevel () const
{
    return m_amr->finestLevel();
}
 
inline 
int 
AmrParGDB::initialBaLevels () const
{
    return Amr::initialBaLevels();
}

inline 
int 
AmrParGDB::subCycle () const
{
    return m_amr->subCycle();
}

inline 
IntVect 
AmrParGDB::refRatio (int level) const
{
    return m_amr->refRatio(level);
}

inline 
int 
AmrParGDB::MaxRefRatio (int level) const
{
    return m_amr->MaxRefRatio(level);
}

#endif
