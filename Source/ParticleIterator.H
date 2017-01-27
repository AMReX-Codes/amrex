#ifndef _PARTICLE_ITERATOR_H_
#define _PARTICLE_ITERATOR_H_

#include <functional>

#include <ParticleContainer.H>

struct PartIterInfo
{
    PartIterInfo (int _level, bool _tiling, const amrex::IntVect& _tilesize)
	: level(_level), tiling(_tiling), tilesize(_tilesize)
	{ }
    PartIterInfo (const PartIterInfo&) = default;
    PartIterInfo (PartIterInfo&&) = default;
    int  level;
    bool tiling;
    amrex::IntVect tilesize;
};

class PartIter
{
public:
    PartIter (WarpXParticleContainer& _mypc, const PartIterInfo& _info);
    PartIter (WarpXParticleContainer& _mypc, PartIterInfo&& _info);
    PartIter (const PartIter&) = delete;
    PartIter (PartIter&&) = delete;
    PartIter& operator=(const PartIter&) = delete;

    bool isValid () { return currentIndex < endIndex; }

    void operator++ () { ++currentIndex; }

    int index () const { return index_map[currentIndex]; }

    int numParticles () const;

    amrex::Box tilebox () const;
    
    amrex::Box validbox () const;

    // apply function to each particle in current tile
    template <typename F>
    void foreach (F&& f)
    {
	auto& pbox = pmap[index_map[currentIndex]];
	if (info.tiling)
	{
	    amrex::Abort("PartIter wip");
	}
	else
	{
	    for (int i = 0, N = pbox.size(); i < N; ++i)
	    {
		f(i,pbox[i]);
	    }
	}
    }

private:
    void Initialize();

    WarpXParticleContainer& mypc;
    WarpXParticleContainer::PMap& pmap;
    PartIterInfo info;

    int currentIndex;
    int beginIndex;
    int endIndex;
    amrex::Array<int> index_map;
};

#endif
