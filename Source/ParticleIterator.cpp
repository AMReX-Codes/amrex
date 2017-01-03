
#include <ParticleIterator.H>

PartIter::PartIter (MyParticleContainer& _mypc, const PartIterInfo& _info)
    : mypc(_mypc),
      pmap(_mypc.GetParticles(_info.level)),
      info(_info)
{
    Initialize();
}

PartIter::PartIter (MyParticleContainer& _mypc, PartIterInfo&& _info)
    : mypc(_mypc),
      pmap(_mypc.GetParticles(_info.level)),
      info(_info)
{
    Initialize();
}

void
PartIter::Initialize ()
{
    if (info.tiling)
    {
	BoxLib::Abort("PartIter wip");
#if 0
	for (auto& pbox : pmap)
	{
	    int np = pbox.size();
	    const IntVect& small = part_box.smallEnd();
	    const IntVect& size  = part_box.size();
	    IntVect ntiles(D_DECL((size[0]+tile_size[0]-1)/tile_size[0],
				  (size[1]+tile_size[1]-1)/tile_size[1],
				  (size[2]+tile_size[2]-1)/tile_size[2]));
	    int ntottile = D_TERM(ntiles[0],*ntiles[1],*ntiles[2]);
	    Array<int> buckets(ntottile);
	    for (int i=0; i<np_box; ++i)
	    {
		const auto& p = orig_pbox[i];
		const auto& cell = p.m_cell;
		IntVect tileindex(D_DECL((cell[0]-small[0])/tile_size[0],
					 (cell[1]-small[1])/tile_size[1],
					 (cell[2]-small[2])/tile_size[2]));
#if (BL_SPACEDIM == 3)
		int idx1d = tileindex[0] + (tileindex[1]-1)*ntiles[0] + (tileindex[2]-1)*ntiles[0]*ntiles[1];
#else
		int idx1d = tileindex[0] + (tileindex[1]-1)*ntiles[0];
#endif
		buckets[idx1d].push_back(i);
	    } 
	}
#endif
    }
    else
    {
	currentIndex = 0;
	beginIndex = 0;
	endIndex = pmap.size();
	for (auto& kv : pmap)
	{
	    index_map.push_back(kv.first);
	}
    }
}

int
PartIter::numParticles () const
{
    const int idx = index_map[currentIndex];
    const auto& pbox = pmap[idx];
    if (info.tiling)
    {
	BoxLib::Abort("PartIter wip");
	return 0;
    }
    else
    {
	return pbox.size();
    }
}

Box
PartIter::tilebox () const
{
    const int idx = index_map[currentIndex];
    const BoxArray& ba = mypc.ParticleBoxArray(info.level);
    return ba[idx]; 
}
    
Box
PartIter::validbox () const
{
    const int idx = index_map[currentIndex];
    const BoxArray& ba = mypc.ParticleBoxArray(info.level);
    return ba[idx]; 
}
