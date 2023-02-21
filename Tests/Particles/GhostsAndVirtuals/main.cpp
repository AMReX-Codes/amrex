#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

#include <array>
#include <iostream>


using namespace amrex;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
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
    Vector<int> rr(nlevs);
    for (int lev = 1; lev < nlevs; lev++)
        rr.at(lev-1) = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    std::array<int, BL_SPACEDIM> is_per;
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per.at(i) = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per.data());
    for (int lev = 1; lev < nlevs; lev++) {
        geom.at(lev).define(amrex::refine(geom.at(lev-1).Domain(), rr.at(lev-1)),
                            &real_box, CoordSys::cartesian, is_per.data());
    }

    Vector<BoxArray> ba(nlevs);
    ba.at(0).define(domain);

    int n_fine = parms.nx;
    for (int lev = 1; lev < nlevs; lev++) {
        n_fine *= rr.at(lev-1);
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba.at(lev).define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba.at(lev).maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap.at(lev) = DistributionMapping{ba.at(lev)};
    }

    typedef AmrParticleContainer<1, 0, 0, 1> MyParticleContainer;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(parms.verbose);

    int num_particles = parms.nppc * AMREX_D_TERM(parms.nx, * parms.ny, * parms.nz);
    bool serialize = true;
    int iseed = 451;
    double mass = 10.0;
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

    int nlevs = 3;

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 64.0);
    }

    int nx,ny,nz;
    nx = ny = nz = 32;
    IntVect domain_lo(AMREX_D_DECL(0 , 0, 0));
    IntVect domain_hi(AMREX_D_DECL(nx - 1, ny - 1, nz-1));
    const Box domain(domain_lo, domain_hi);

    amrex::Print()<<"Ascii test always uses nx=ny=nz=32, nlevs=3, ProbHi=64"<<std::endl;

    // Define the refinement ratio
    Vector<int> rr(nlevs);
    rr[0] = 2;
    for (int lev = 1; lev < nlevs; lev++)
        rr.at(lev) = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    std::array<int, BL_SPACEDIM> is_per;
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per.data());
    for (int lev = 1; lev < nlevs; lev++) {
        geom.at(lev).define(amrex::refine(geom.at(lev-1).Domain(), rr.at(lev-1)),
                            &real_box, CoordSys::cartesian, is_per.data());
    }

    const std::string regrid_grids_file = "fixed_grids.init";

    Vector<BoxArray> ba(nlevs);
    ba.at(0).define(domain);

    //Create BoxArray similar to InitAmr in AMReX_Amr.cpp
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
                bx.refine(rr.at(lev-1));

                bl.push_back(bx);
            }
            ba.at(lev-1).define(bl);
        }
        is.close();
#undef STRIP
    }

    //Create "original" BoxArray based on a hierarchical grid
    Vector<BoxArray> ba_orig(nlevs);
    ba_orig.at(0).define(domain);

    int n_fine = parms.nx;
    for (int lev = 1; lev < nlevs; lev++) {
        n_fine *= rr.at(lev-1);
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba_orig.at(lev).define(refined_patch);
    }

    // break the BoxArrays at all levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba.at(lev).maxSize(parms.max_grid_size);
        ba_orig.at(lev).maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    Vector<DistributionMapping> dmap_orig(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap.at(lev) = DistributionMapping{ba.at(lev)};
        dmap_orig.at(lev) = DistributionMapping{ba_orig.at(lev)};
    }

    typedef AmrParticleContainer<4, 0, 0, 0> MyParticleContainer;
    using PType = typename MyParticleContainer::SuperParticleType;
    MyParticleContainer myPC(geom, dmap_orig, ba_orig, rr);

    myPC.SetVerbose(parms.verbose);

    //Initialize particles from ascii file with extradata=4
    myPC.InitFromAsciiFile("particle_file.init", 4);

    //Regrid to the more compilcated BoxArray similar to NeighborParticleContainer Regrid
    for (int lev = 0; lev < nlevs; lev++) {
        myPC.SetParticleBoxArray(lev, ba.at(lev));
        myPC.SetParticleDistributionMap(lev, dmap.at(lev));
    }
    myPC.Redistribute();

    //Check created particle containers with a summation and number of particles
    Real tol = 1e-4;
    {

        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT(virtPC.TotalNumberOfParticles(true,false)==3);
        AMREX_ALWAYS_ASSERT(std::abs(sum_test-3029.00000028022578)<tol);

    }

    {

        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(2, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT(virtPC.TotalNumberOfParticles(true,false)==0);
        AMREX_ALWAYS_ASSERT(std::abs(sum_test-0.0)<tol);
    }

    {
        const int ngrow = 0;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT(ghostPC.TotalNumberOfParticles(true,false)==0);
        AMREX_ALWAYS_ASSERT(std::abs(sum_test-0.0)<tol);
    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT(ghostPC.TotalNumberOfParticles(true,false)==3);
        AMREX_ALWAYS_ASSERT(std::abs(sum_test-3035.00000001795206)<tol);
    }

    {
        const int ngrow = 4;
        const int src_lev = 1;
        const int dst_lev = 2;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT(ghostPC.TotalNumberOfParticles(true,false)==1);
        AMREX_ALWAYS_ASSERT(std::abs(sum_test-1005.00000009692667)<tol);
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
    Vector<int> rr(nlevs);
    for (int lev = 1; lev < nlevs; lev++)
        rr.at(lev-1) = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    std::array<int, BL_SPACEDIM> is_per;
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per.at(i) = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom.at(0).define(domain, &real_box, CoordSys::cartesian, is_per.data());
    for (int lev = 1; lev < nlevs; lev++) {
        geom.at(lev).define(amrex::refine(geom.at(lev-1).Domain(), rr.at(lev-1)),
                            &real_box, CoordSys::cartesian, is_per.data());
    }

    Vector<BoxArray> ba(nlevs);
    ba.at(0).define(domain);

    int n_fine = parms.nx;
    for (int lev = 1; lev < nlevs; lev++) {
        n_fine *= rr.at(lev-1);
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba.at(lev).define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba.at(lev).maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap.at(lev) = DistributionMapping{ba.at(lev)};
    }

    typedef AmrParticleContainer<4, 0, 0, 0> MyParticleContainer;
    using PType = typename MyParticleContainer::SuperParticleType;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc;
    int iseed = 451;
    double mass = 10.0;
    double xvel, yvel, zvel;
    Real total_virts_test = 0.0;
    Real tol = 1e-4;
    xvel = 1.0;
    yvel = 2.0;
    zvel = 3.0;
    total_virts_test+=mass+xvel+yvel+zvel;
    MyParticleContainer::ParticleInitData pdata = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitRandomPerBox(num_particles, iseed, pdata);

    mass = 1000.0;
    xvel = 1.0e-5;
    yvel = 2.0e-5;
    zvel = 3.0e-5;
    total_virts_test+=mass+xvel+yvel+zvel;
    MyParticleContainer::ParticleInitData pdata_big = {{mass, xvel, yvel, zvel},{}, {}, {}};
    MyParticleContainer myPC_tmp(geom, dmap, ba, rr);
    myPC_tmp.InitRandomPerBox(num_particles, iseed, pdata_big);
    myPC.addParticles(myPC_tmp);
    myPC_tmp.clearParticles();

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of virts: "<<sum_test<<std::endl;
        AMREX_ALWAYS_ASSERT((virtPC.AggregationType()=="None" && std::abs(sum_test - total_virts_test * parms.nx * parms.ny * parms.nz / (32 * 32 * 32) * (16 * 16 * 16) / (parms.max_grid_size * parms.max_grid_size * parms.max_grid_size) * parms.nppc * parms.nppc * parms.nppc) < tol) || ParallelDescriptor::NProcs() % 2 != 0 || virtPC.AggregationType() == "Cell");

    }

    {
        const int ngrow = 1;
        const int src_lev = 0;
        const int dst_lev = 1;
        MyParticleContainer ghostPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType ghosts;
        myPC.CreateGhostParticles(src_lev, ngrow, ghosts);
        ghostPC.AddParticlesAtLevel(ghosts, dst_lev, ngrow);
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
    }

    mass = 1000000.0;
    xvel = 1.0e-8;
    yvel = 2.0e-8;
    zvel = 3.0e-8;
    MyParticleContainer::ParticleInitData pdata_bigger = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC_tmp.InitRandomPerBox(num_particles, iseed+5, pdata_bigger);
    myPC.addParticles(myPC_tmp);
    myPC_tmp.clearParticles();

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
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
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        Long id_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Long
            {
                return p.id();
            }
        );
        amrex::ParallelAllReduce::Sum(id_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
        amrex::Print().SetPrecision(18)<<"Found sum of id of ghosts: "<<id_test<<" ?= "<<ghostPC.TotalNumberOfParticles(true,false)*GhostParticleID<<std::endl;
        AMREX_ALWAYS_ASSERT(id_test==ghostPC.TotalNumberOfParticles(true,false)*GhostParticleID);
        AMREX_ALWAYS_ASSERT(ghostPC.TotalNumberOfParticles(true,false)==ghostPC.TotalNumberOfParticles(false,false));
    }
}

void test_ghosts_and_virtuals_onepercell (TestParams& parms)
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
    Vector<int> rr(nlevs);
    for (int lev = 1; lev < nlevs; lev++)
        rr.at(lev-1) = 2;

    // This sets the boundary conditions to be doubly or triply periodic
    std::array<int, BL_SPACEDIM> is_per;
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per.at(i) = 1;

    // This defines a Geometry object which is useful for writing the plotfiles
    Vector<Geometry> geom(nlevs);
    geom.at(0).define(domain, &real_box, CoordSys::cartesian, is_per.data());
    for (int lev = 1; lev < nlevs; lev++) {
        geom.at(lev).define(amrex::refine(geom.at(lev-1).Domain(), rr.at(lev-1)),
                            &real_box, CoordSys::cartesian, is_per.data());
    }

    Vector<BoxArray> ba(nlevs);
    ba.at(0).define(domain);

    int n_fine = parms.nx;
    for (int lev = 1; lev < nlevs; lev++) {
        n_fine *= rr.at(lev-1);
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba.at(lev).define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba.at(lev).maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap.at(lev) = DistributionMapping{ba.at(lev)};
    }

    typedef AmrParticleContainer<4, 0, 0, 0> MyParticleContainer;
    using PType = typename MyParticleContainer::SuperParticleType;
    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    double mass = 10.0;
    double xvel, yvel, zvel;
    Real xoff, yoff, zoff;
    xvel = 1.0;
    yvel = 2.0;
    zvel = 3.0;
    xoff = 0.5;
    yoff = 0.5;
    zoff = 0.5;
    MyParticleContainer::ParticleInitData pdata = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitOnePerCell(xoff, yoff, zoff, pdata);

    mass = 1000.0;
    xvel = 1.0e-5;
    yvel = 2.0e-5;
    zvel = 3.0e-5;
    xoff = 0.25;
    yoff = 0.25;
    zoff = 0.25;
    MyParticleContainer::ParticleInitData pdata_big = {{mass, xvel, yvel, zvel},{}, {}, {}};
    MyParticleContainer myPC_tmp(geom, dmap, ba, rr);
    myPC_tmp.InitOnePerCell(xoff, yoff, zoff, pdata_big);
    myPC.addParticles(myPC_tmp);
    myPC_tmp.clearParticles();

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
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
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
        amrex::Print().SetPrecision(18)<<"Found sum of ghosts: "<<sum_test<<std::endl;
    }

    mass = 1000000.0;
    xvel = 1.0e-8;
    yvel = 2.0e-8;
    zvel = 3.0e-8;
    xoff = 0.75;
    yoff = 0.75;
    zoff = 0.75;
    MyParticleContainer::ParticleInitData pdata_bigger = {{mass, xvel, yvel, zvel},{}, {}, {}};
    myPC.InitOnePerCell(xoff, yoff, zoff, pdata_bigger);

    {
        MyParticleContainer virtPC(geom, dmap, ba, rr);
        MyParticleContainer::ParticleTileType virts;
        myPC.CreateVirtualParticles(1, virts);
        virtPC.AddParticlesAtLevel(virts, 0);
        Real sum_test = amrex::ReduceSum(virtPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
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
        Real sum_test = amrex::ReduceSum(ghostPC,
            [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
            {
                return std::abs(p.rdata(0)) +
                       std::abs(p.rdata(1)) +
                       std::abs(p.rdata(2)) +
                       std::abs(p.rdata(3));
            }
        );
        amrex::ParallelAllReduce::Sum(sum_test,ParallelContext::CommunicatorSub());
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

  amrex::Print()<<"Ascii test"<<std::endl;

  test_ghosts_and_virtuals_ascii(parms);

  amrex::Print()<<"Original test"<<std::endl;
  test_ghosts_and_virtuals(parms);

#ifndef AMREX_USE_SYCL
  amrex::Print()<<"RandomPerBox test"<<std::endl;
  test_ghosts_and_virtuals_randomperbox(parms);
#endif

  amrex::Print()<<"OnePerCell test"<<std::endl;
  test_ghosts_and_virtuals_onepercell(parms);

  amrex::Finalize();
}
