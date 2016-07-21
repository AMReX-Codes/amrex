
#include <MultiMask.H>

MultiMask::MultiMask () { }

MultiMask::MultiMask (const BoxArray&            ba,
		      const DistributionMapping& dm,
		      int                        ncomp)
    : FabArray<Mask>(ba,ncomp,0,dm)
{
}


