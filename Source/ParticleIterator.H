#ifndef _PARTICLE_ITERATOR_H_
#define _PARTICLE_ITERATOR_H_

#include <functional>

#include <ParticleContainer.H>

struct PartIterInfo
{
    PartIterInfo (int _level, bool _tiling, const IntVect& _tilesize)
	: level(_level), tiling(_tiling), tilesize(_tilesize)
	{ }
    PartIterInfo (const PartIterInfo&) = default;
    PartIterInfo (PartIterInfo&&) = default;
    int  level;
    bool tiling;
    IntVect tilesize;
};

class PartIter
{
public:
    PartIter (MyParticleContainer& _mypc, const PartIterInfo& _info);
    PartIter (MyParticleContainer& _mypc, PartIterInfo&& _info);
    PartIter (const PartIter&) = delete;
    PartIter (PartIter&&) = delete;
    PartIter& operator=(const PartIter&) = delete;

    bool isValid () { return currentIndex < endIndex; }

    void operator++ () { ++currentIndex; }

    int index () const { return index_map[currentIndex]; }

    int numParticles () const;

    Box tilebox () const;
    
    Box validbox () const;

    // apply function to each particle in current tile
    void foreach (std::function<void(int, MyParticleContainer::ParticleType&)> f);

private:
    void Initialize();

    MyParticleContainer& mypc;
    MyParticleContainer::PMap& pmap;
    PartIterInfo info;

    int currentIndex;
    int beginIndex;
    int endIndex;
    Array<int> index_map;
};

#endif
