#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include "AMReX_Particles.H"
#include "AMReX_PlotFileUtil.H"

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  bool verbose;
};

void test_init_ascii (TestParams& parms)
{
    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    
    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);
    
    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++) 
        is_per[i] = 1;
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(parms.max_grid_size);

    DistributionMapping dmap(ba);

    typedef ParticleContainer<1, 0, AMREX_SPACEDIM> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba);

    myPC.InitFromAsciiFile("particles.txt", 1 + AMREX_SPACEDIM);

    // should be 8
    amrex::Print() << myPC.TotalNumberOfParticles() << "\n";

    // should be 8010.0
    using PType = MyParticleContainer::SuperParticleType;
    amrex::Print() << amrex::ReduceSum(myPC,
                                       [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                       {
                                           Real total = 0.0;
                                           for (int i = 0; i < 1 + AMREX_SPACEDIM; ++i)
                                           {
                                               total += p.rdata(i);
                                           }
                                           return total;
                                       })
                   << "\n";
}

int main(int argc, char* argv[])
{
  amrex::Initialize(argc,argv);
  
  ParmParse pp;
  
  TestParams parms; 
  pp.get("nx", parms.nx);
  pp.get("ny", parms.ny);
  pp.get("nz", parms.nz);
  pp.get("max_grid_size", parms.max_grid_size);
  
  test_init_ascii(parms);
  
  amrex::Finalize();
}
