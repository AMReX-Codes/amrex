
#include <ParticleContainer.H>
#include <MultiFab.H>

MyParticleContainer::MyParticleContainer (const Array<Geometry>            & geom,
					  const Array<DistributionMapping> & dmap,
					  const Array<BoxArray>            & ba,
					  const Array<int>                 & rr)
    : ParticleContainer<PIdx::nattribs,0> (geom,dmap,ba,rr)
{
    MultiFab foo(ba[0],1,0,dmap[0],Fab_noallocate);
    for (MFIter mfi(foo); mfi.isValid(); ++mfi)
    {
	int i = mfi.index();
	partdata[i] = Array<std::unique_ptr<Array<Real> > > (PIdx::npartdata);
	for (auto& d : partdata[i])
	{
	    d = std::unique_ptr<Array<Real> >(new Array<Real>());
	}
    }
}
