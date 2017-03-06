#include <iostream>
#include <map>
#include <vector>
#include <type_traits>

#include "AMReX_Array.H"
#include "AMReX_Particles.H"

using namespace amrex;

void test_POD()
{
  const int NR = 4;
  const int NI = 2;

  std::cout << "Is POD: "      << std::is_pod<Particle<NR, NI> >::value << std::endl;
  std::cout << "Is trivial: "  << std::is_trivial<Particle<NR, NI> >::value << std::endl;
  std::cout << "Is standard: " << std::is_standard_layout<Particle<NR, NI> >::value << std::endl;
  std::cout << "Size: "        << sizeof(Particle<NR, NI>) << std::endl;

  Particle<NR, NI> p;
  p.id() = 1;
  p.cpu() = 0;
  p.pos(0) = 1.0;
  p.pos(1) = 2.0;
  p.pos(2) = 3.0;
  p.rdata(0) = 4.0;
  p.rdata(1) = -1.0;
  p.rdata(2) = -2.0;
  p.rdata(3) = -3.0;

  p.idata(0) = -7;
  p.idata(1) =  9;

  std::cout << p.id() << " " << p.cpu() << std::endl;
  std::cout << p.pos(0) << " " << p.pos(1) << " " << p.pos(2) << std::endl;
  std::cout << p.rdata(0) << " " << p.rdata(1) << " " << p.rdata(2) << " " << p.rdata(3) << std::endl;
  std::cout << p.idata(0) << " " << p.idata(1) << std::endl;
}

int main(int argc, char* argv[])
{

  amrex::Initialize(argc,argv);

  test_POD();

  int ncell = 64;
  int nlevs = 2;
  int coord = 0;

  RealBox real_box, fine_box;
  for (int n = 0; n < BL_SPACEDIM; n++)
    {
      real_box.setLo(n,0.0);
      real_box.setHi(n,1.0);
      fine_box.setLo(n,0.4);
      fine_box.setHi(n,0.6);
    }

  IntVect domain_lo(0 , 0, 0); 
  IntVect domain_hi(ncell-1, ncell-1, ncell-1); 
 
  const Box domain(domain_lo, domain_hi);

  Array<int> rr(nlevs-1);
  for (int lev = 1; lev < nlevs; lev++)
    rr[lev-1] = 2;
 
  int is_per[BL_SPACEDIM];
  for (int i = 0; i < BL_SPACEDIM; i++) is_per[i] = 1;

  Array<Geometry> geom(nlevs);
  geom[0].define(domain, &real_box, coord, is_per);
  geom[1].define(amrex::refine(geom[0].Domain(), rr[0]),
		 &real_box, coord, is_per);

  Array<BoxArray> ba(nlevs);
  ba[0].define(domain);  

  int nfine = ncell*rr[0];
  IntVect refined_lo(nfine/4,nfine/4,nfine/4); 
  IntVect refined_hi(3*nfine/4-1,3*nfine/4-1,3*nfine/4-1);

  // Build a box for the level 1 domain
  Box refined_patch(refined_lo, refined_hi);
  ba[1].define(refined_patch);

  int max_grid_size = 32;
  for (int lev = 0; lev < nlevs; lev++)
    ba[lev].maxSize(max_grid_size);

  Array<DistributionMapping> dmap(nlevs);
  for (int lev = 0; lev < nlevs; lev++)
    dmap[lev].define(ba[lev]);

  typedef ParticleContainer<1+BL_SPACEDIM, 0, 5> MyParticleContainer;

  // Build a new particle container to hold my particles.
  std::unique_ptr<MyParticleContainer> MyPC(new MyParticleContainer(geom, dmap, ba, rr));
  MyPC->communicate_comp[2] = false;
  MyPC->communicate_comp[3] = false;
  MyPC->communicate_comp[4] = false;

  MFInfo Fab_noallocate;
  Fab_noallocate.SetAlloc(false);
  MultiFab dummy_mf(ba[0], dmap[0], 1, 0, Fab_noallocate);

  MyPC->InitOnePerCell(0.5, 0.5, 0.5, 1.0, dummy_mf);

  for (int i = 0; i < 10; i ++) {
    MyPC->MoveRandom();
  }

  MyPC->WriteAsciiFile("particles");
    
  amrex::Finalize();
}
