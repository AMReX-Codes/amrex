#ifndef _ParGDB_H
#define _ParGDB_H

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

    virtual const BoxArray& boxArray (int level) const = 0;
    virtual const MultiFab& multiFab (int level) const = 0;

    virtual bool LevelDefined (int level) const = 0;
    virtual int finestLevel () const = 0;
    virtual int initialBaLevels () const = 0;

    virtual int subCycle () const = 0;
    virtual IntVect refRatio (int level) const = 0;
    virtual int MaxRefRatio (int level) const = 0;

    virtual bool ParticlesOnSameGrids (int level) const = 0;
};

#if 0
// we use this for non-Amr particle code
class ParGDB
    : public ParGDBBase
{
public:

    virtual const Geometry& Geom (int level) const;
    virtual const DistributionMapping& ParticleDistributionMap 
                                             (int level) const;
    virtual const BoxArray& ParticleBoxArray (int level) const;

    virtual const BoxArray& boxArray (int level) const;
    virtual const MultiFab& multiFab (int level) const;

    virtual bool LevelDefined (int level) const;
    virtual int finestLevel () const;
    virtual int initialBaLevels () const;

    virtual int subCycle () const;
    virtual IntVect refRatio (int level) const;
    virtual int MaxRefRatio (int level) const;

    virtual bool ParticlesOnSameGrids (int level) const;

protected:

    // we are going to put some date here, like geom, dm and ba

};
#endif

#endif
