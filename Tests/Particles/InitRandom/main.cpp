#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Particles.H>

#include <cstdio>

using namespace amrex;

void set_grids_nested (Vector<Box>& domains,
                       Vector<BoxArray>& grids,
                       Vector<IntVect>& ref_ratio);
void test ();
void testSOA ();

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    test();
    testSOA();
    amrex::Finalize();
}

void testSOA ()
{
    int ncells, max_grid_size, nlevs, nppc;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);

    Vector<Box> domains;
    Vector<BoxArray> ba;
    Vector<IntVect> ref_ratio;

    set_grids_nested(domains, ba, ref_ratio);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[] = {AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object for each level
    Vector<Geometry> geom(nlevs);
    geom[0].define(domains[0], &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(domains[lev], &real_box, CoordSys::cartesian, is_per);
    }

    Vector<DistributionMapping> dmap(nlevs);

    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }

    // Add some particles
    constexpr int NReal  = 12;
    constexpr int NInt   = 4;

    using MyPC = ParticleContainerPureSoA<NReal, NInt>;
    MyPC myPC(geom, dmap, ba, ref_ratio);
    myPC.SetVerbose(false);

    int num_particles = nppc * AMREX_D_TERM(ncells, * ncells, * ncells);
    bool serialize = false;
    int iseed = 451;
    MyPC::ParticleInitData pdata = {{}, {},
                      {1.0, 2.0, 3.0, 4.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0},
                      {5, 14, 15, 16}};

    myPC.InitRandom(num_particles, iseed, pdata, serialize);
    amrex::Print() << "Generated " << myPC.TotalNumberOfParticles() << " particles. \n";
}

void test ()
{
    int ncells, max_grid_size, nlevs, nppc;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);

    Vector<Box> domains;
    Vector<BoxArray> ba;
    Vector<IntVect> ref_ratio;

    set_grids_nested(domains, ba, ref_ratio);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[] = {AMREX_D_DECL(1,1,1)};

    // This defines a Geometry object for each level
    Vector<Geometry> geom(nlevs);
    geom[0].define(domains[0], &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(domains[lev], &real_box, CoordSys::cartesian, is_per);
    }

    Vector<DistributionMapping> dmap(nlevs);

    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
    }

    // Add some particles
    constexpr int NStructReal = 4;
    constexpr int NStructInt  = 1;
    constexpr int NArrayReal  = 8;
    constexpr int NArrayInt   = 3;

    using MyPC = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    MyPC myPC(geom, dmap, ba, ref_ratio);
    myPC.SetVerbose(false);

    int num_particles = nppc * AMREX_D_TERM(ncells, * ncells, * ncells);
    bool serialize = false;
    int iseed = 451;
    MyPC::ParticleInitData pdata = {{1.0, 2.0, 3.0, 4.0},
                                    {5},
                                    {6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0},
                                    {14, 15, 16}};

    myPC.InitRandom(num_particles, iseed, pdata, serialize);
    amrex::Print() << "Generated " << myPC.TotalNumberOfParticles() << " particles. \n";
}

void set_grids_nested (Vector<Box>& domains,
                       Vector<BoxArray>& grids,
                       Vector<IntVect>& ref_ratio)
{
    int ncells, max_grid_size, nlevs;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("nlevs", nlevs);

    AMREX_ALWAYS_ASSERT(nlevs < 2); // relax this later

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(ncells-1, ncells-1, ncells-1));

    domains.resize(nlevs);
    domains[0].setSmall(domain_lo);
    domains[0].setBig(domain_hi);

    ref_ratio.resize(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++) {
        ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));
    }

    grids.resize(nlevs);
    grids[0].define(domains[0]);

    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = ncells*ref_ratio[0][0];
        IntVect refined_lo(AMREX_D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(AMREX_D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        grids[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        grids[lev].maxSize(max_grid_size);
    }

    for (int lev = 1; lev < nlevs; lev++) {
        domains[lev] = amrex::refine(domains[lev-1], ref_ratio[lev-1]);
    }
}
