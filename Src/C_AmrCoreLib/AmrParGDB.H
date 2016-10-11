#ifndef _AmrParGDB_H
#define _AmrParGDB_H

#include <ParGDB.H>
#include <AmrCore.H>

class AmrParGDB
    : public ParGDBBase
{
public:

    AmrParGDB (AmrCore* amr)
	: m_amrcore(amr) { BL_ASSERT(amr != 0); }

    virtual ~AmrParGDB () {;}
    
    virtual const Geometry& Geom (int level) const override;
    virtual const DistributionMapping& ParticleDistributionMap 
                                             (int level) const override;
    virtual const BoxArray& ParticleBoxArray (int level) const override;
    virtual const BoxArray&         boxArray (int level) const override;

    virtual void SetParticleBoxArray (int level,
				      const DistributionMapping& new_dm,
				      const BoxArray           & new_ba) override;

    virtual bool LevelDefined (int level) const override;
    virtual int finestLevel () const override;
    virtual int initialBaLevels () const override;

    virtual int subCycle () const override;
    virtual IntVect refRatio (int level) const override;
    virtual int MaxRefRatio (int level) const override;

protected:
    AmrCore* m_amrcore;
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
    return m_amrcore->ParticleDistributionMap(level);
}

inline 
const BoxArray&
AmrParGDB::ParticleBoxArray (int level) const
{
    return m_amrcore->ParticleBoxArray(level);
}

inline
const BoxArray&
AmrParGDB::boxArray (int level) const
{
    return m_amrcore->boxArray(level);
}

inline
void AmrParGDB::SetParticleBoxArray (int level,
				     const DistributionMapping& new_dm,
				     const BoxArray           & new_ba)
{
    // xxxxx
    BoxLib::Abort("AmrParGDB::SetParticleBoxArray: not implemented yet");
}

inline 
bool 
AmrParGDB::LevelDefined (int level) const
{
//    return m_amrcore->getAmrLevels().defined(level);
    // xxxxx
    return true;
}

inline 
int 
AmrParGDB::finestLevel () const
{
    return m_amrcore->finestLevel();
}
 
inline 
int 
AmrParGDB::initialBaLevels () const
{
    // xxxxx
//    return Amr::initialBaLevels();
    return 0;
}

inline 
int 
AmrParGDB::subCycle () const
{
    // xxxxx
//    return m_amrcore->subCycle();
    return 0;
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

#endif
