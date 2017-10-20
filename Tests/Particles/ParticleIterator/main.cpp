#include <iostream>
#include <map>
#include <vector>

#include <AMReX_Vector.H>
#include "AMReX_FabArray.H"
#include "AMReX_Particles.H"

using namespace amrex;

int main(int argc, char* argv[])
{

  amrex::Initialize(argc,argv);

  int ncell = 48;
  int max_grid_size = 32;
  int nlevs = 1;
  int coord = 0;

  RealBox real_box;
  for (int n = 0; n < BL_SPACEDIM; n++)
    {
      real_box.setLo(n,0.0);
      real_box.setHi(n,1.0);
    }

  IntVect domain_lo(0 , 0, 0); 
  IntVect domain_hi(ncell-1, ncell-1, ncell-1); 
 
  const Box domain(domain_lo, domain_hi);

  Vector<int> rr(nlevs-1);
  for (int lev = 1; lev < nlevs; lev++)
    rr[lev-1] = 2;
 
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1;

  Vector<Geometry> geom(nlevs);
  geom[0].define(domain, &real_box, coord, is_per);

  Vector<BoxArray> ba(nlevs);
  ba[0].define(domain);  

  for (int lev = 0; lev < nlevs; lev++)
    ba[lev].maxSize(max_grid_size);

  Vector<DistributionMapping> dmap(nlevs);
  for (int lev = 0; lev < nlevs; lev++)
    dmap[lev].define(ba[lev]);

  typedef ParticleContainer<1+BL_SPACEDIM> MyParticleContainer;
  MyParticleContainer MyPC(geom, dmap, ba, rr);

  MyParticleContainer::ParticleInitData pdata = {1.0};
  MyPC.InitOnePerCell(0.5, 0.5, 0.5, pdata);
  MyPC.do_tiling = true;

  amrex::AllPrintToFile("outside") << "outside parallel region. \n";

#ifdef _OPENMP
#pragma omp parallel
#endif
  for (ParIter<1+BL_SPACEDIM> mfi(MyPC, 0); mfi.isValid(); ++mfi) {
      amrex::AllPrintToFile("particle_iterator_out") << mfi.index() << " " << mfi.tileIndex() << "\n";
  }
      
  amrex::Finalize();
}
