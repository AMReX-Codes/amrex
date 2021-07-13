#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
  std::string regrid_file;
};

void test_ghosts_and_virtuals (TestParams& parms)
{
    int nlevs = parms.nlevs;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }

    IntVect domain_lo(AMREX_D_DECL(0 , 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);

    for (int lev = 1; lev < nlevs; lev++) {
        int n_fine = parms.nx*rr[lev-1];
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[lev].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }

    typedef AmrParticleContainer<1, 0, 0, 0> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(parms.verbose);

    int num_particles = parms.nppc * AMREX_D_TERM(parms.nx, * parms.ny, * parms.nz);
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    MyParticleContainer::ParticleInitData pdata = {{mass}, {}, {}, {}};

    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
    }
}

void test_ghosts_and_virtuals_ascii (TestParams& parms)
{

    int nlevs = parms.nlevs;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 64.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0 , 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    rr[0] = 2;
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev] = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    const auto regrid_grids_file = parms.regrid_file;

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    if (nlevs > 0 && !regrid_grids_file.empty())
    {
#define STRIP while( is.get() != '\n' ) {}
        std::ifstream is(regrid_grids_file.c_str(),std::ios::in);

        if (!is.good())
            amrex::FileOpenFailed(regrid_grids_file);

        int in_finest,ngrid;

        is >> in_finest;
        STRIP;
        ba.resize(in_finest);
        for (int lev = 2; lev <= in_finest; lev++)
        {
            BoxList bl;
            is >> ngrid;
            STRIP;
            for (int i = 0; i < ngrid; i++)
            {
                Box bx;
                is >> bx;
                STRIP;
                 bx.refine(rr[lev-1]);
         /*
                 for (int idim = 0 ; idim < AMREX_SPACEDIM; ++idim)
                 {
                     if (bx.length(idim) > parms.max_grid_size)
                     {
                         std::ostringstream ss;
                         ss << "Grid " << bx << " too large" << '\n';
                         amrex::Error(ss.str());
                     }
                 }*/
                 bl.push_back(bx);
            }
            ba[lev-1].define(bl);
        //            amrex::Print()<<(geom[lev-1].ProbHi(0))<<"\n"<<ba[lev-1]<<std::endl;
        }
        is.close();
#undef STRIP
    }
    /*
    for (int lev = 0; lev < nlevs; lev++)
    AmrMesh::ChopGrids(lev,ba[lev],ParallelDescriptor::NProcs());
    */
    Vector<BoxArray> ba_orig(nlevs);
    ba_orig[0].define(domain);

    int n_fine = parms.nx;
    for (int lev = 1; lev < nlevs; lev++) {
        n_fine *= rr[lev-1];
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba_orig[lev].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
        ba_orig[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    Vector<DistributionMapping> dmap_orig(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        dmap_orig[lev] = DistributionMapping{ba_orig[lev]};
    }

    typedef AmrParticleContainer<4, 0, 0, 0> MyParticleContainer;
    using PType = typename amrex::ParticleContainer<4, 0>::SuperParticleType;
    MyParticleContainer myPC(geom, dmap_orig, ba_orig, rr);

    myPC.SetVerbose(parms.verbose);
    /*
    for (int lev = 0; lev < nlevs; lev++) {
    amrex::Print()<<myPC.ParticleBoxArray(lev)<<std::endl;
    amrex::Print()<<myPC.ParticleDistributionMap(lev)<<std::endl;
    }*/

    myPC.InitFromAsciiFile("particle_file.init", 4, nullptr);

    for (int lev = 0; lev < nlevs; lev++) {
    myPC.SetParticleBoxArray(lev, ba[lev]);
    myPC.SetParticleDistributionMap(lev, dmap[lev]);
    }
    myPC.Redistribute();

    {

        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
        virtPC.WriteAsciiFile("virt_1_0");
    }

    {

        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(2, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
        virtPC.WriteAsciiFile("virt_2_0");
    }

    {
        const int ngrow = 0;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        ghostPC.WriteAsciiFile("ghost_0_0_1");
    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        ghostPC.WriteAsciiFile("ghost_1_0_1");
    }

    {
        const int ngrow = 4;
        const int src_lev = 1;
        const int dst_lev = 2;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        ghostPC.WriteAsciiFile("ghost_4_1_2");
    }
}

void test_ghosts_and_virtuals_randomperbox (TestParams& parms)
{

    int nlevs = parms.nlevs;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    RealBox fine_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
       fine_box.setLo(n,0.25);
       fine_box.setHi(n,0.75);
    }

    IntVect domain_lo(AMREX_D_DECL(0 , 0, 0));
    IntVect domain_hi(AMREX_D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<int> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);

    if (nlevs > 1) {
        int n_fine = parms.nx*rr[0];
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }

    typedef AmrParticleContainer<4, 0, 0, 0> MyParticleContainer;
    using PType = typename amrex::ParticleContainer<4, 0>::SuperParticleType;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc;
    int iseed = 451;
    Real mass = 10.0;
    Real xvel, yvel, zvel;
    xvel = 1.0;
    yvel = 2.0;
    zvel = 3.0;
    MyParticleContainer::ParticleInitData pdata = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitRandomPerBox(num_particles, iseed, pdata);

    mass = 1000.0;
    xvel = 1.0e-5;
    yvel = 2.0e-5;
    zvel = 3.0e-5;
    MyParticleContainer::ParticleInitData pdata_big = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitRandomPerBox(num_particles, iseed, pdata_big);

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;

    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
    amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
    }

    mass = 1000000.0;
    xvel = 1.0e-8;
    yvel = 2.0e-8;
    zvel = 3.0e-8;
    MyParticleContainer::ParticleInitData pdata_bigger = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitRandomPerBox(num_particles, iseed+5, pdata_bigger);

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real { return (amrex::Math::abs(p.rdata(0))+amrex::Math::abs(p.rdata(1))+amrex::Math::abs(p.rdata(2))+amrex::Math::abs(p.rdata(3))); });
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
    }
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
  pp.get("nlevs", parms.nlevs);
  pp.get("nppc", parms.nppc);
  if (parms.nppc < 1 && ParallelDescriptor::IOProcessor())
    amrex::Abort("Must specify at least one particle per cell");

  parms.verbose = false;
  pp.query("verbose", parms.verbose);

  if (parms.verbose && ParallelDescriptor::IOProcessor()) {
    std::cout << std::endl;
    std::cout << "Number of particles per cell : ";
    std::cout << parms.nppc  << std::endl;
    std::cout << "Size of domain               : ";
    std::cout << "Num levels: ";
    std::cout << parms.nlevs << std::endl;
    std::cout << parms.nx << " " << parms.ny << " " << parms.nz << std::endl;
  }

  parms.regrid_file = "";
  pp.query("regrid_file", parms.regrid_file);
  std::cout<<parms.regrid_file<<std::endl;
  test_ghosts_and_virtuals_ascii(parms);

  //  test_ghosts_and_virtuals(parms);

  amrex::Finalize();
}
