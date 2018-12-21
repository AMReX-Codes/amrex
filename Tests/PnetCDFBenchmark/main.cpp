#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Particles.H>

#ifdef BL_PNETCDF
#include <WritePlotfilePnetCDF.H>
#endif

using namespace amrex;

int main(int argc, char* argv[])
{    
    amrex::Initialize(argc,argv);
    {
    const int nghost = 0;
    int ncells, max_grid_size, ncomp, nlevs, nppc;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("ncomp", ncomp);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);
    
    AMREX_ALWAYS_ASSERT(nlevs < 2); // relax this later

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(ncells-1, ncells-1, ncells-1)); 
    const Box domain(domain_lo, domain_hi);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // Define the refinement ratio
    Vector<IntVect> ref_ratio(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object for each level
    Vector<Geometry> geom(nlevs);
    geom[0].define(domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
	geom[lev].define(amrex::refine(geom[lev-1].Domain(), ref_ratio[lev-1]),
			 &real_box, CoordSys::cartesian, is_per);
    }
    
    Vector<BoxArray> ba(nlevs);
    ba[0].define(domain);
    
    // Now we make the refined level be the center eighth of the domain
    if (nlevs > 1) {
        int n_fine = ncells*ref_ratio[0][0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4)); 
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }
    
    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);

    Vector<std::unique_ptr<MultiFab> > mf(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        mf[lev].reset(new MultiFab(ba[lev], dmap[lev], ncomp, nghost));
        mf[lev]->setVal(lev);
    }

    // Add some particles
    constexpr int NStructReal = 4;
    constexpr int NStructInt  = 1;
    constexpr int NArrayReal  = 8;
    constexpr int NArrayInt   = 3;

    typedef ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt> MyPC;
    MyPC myPC(geom, dmap, ba, ref_ratio);
    myPC.SetVerbose(false);

    int num_particles = nppc * AMREX_D_TERM(ncells, * ncells, * ncells);
    bool serialize = false;
    int iseed = 451;
    MyPC::ParticleInitData pdata = {1.0, 2.0, 3.0, 4.0, 5, 6.0,
                                    7.0, 8.0, 9.0, 10.0, 11.0,
                                    12.0, 13.0, 14, 15, 16};
    
    
    myPC.InitRandom(num_particles, iseed, pdata, serialize);
    
    // these don't really matter, make something up
    const Real time = 0.0;
    const Real dt = 0.0;

    Vector<std::string> varnames;
    for (int i = 0; i < ncomp; ++i)
    {
        varnames.push_back("component_" + std::to_string(i));
    }

#ifdef BL_PNETCDF   
    WriteMultiLevelPlotfilePNETCDF("plt00000", nlevs, amrex::GetVecOfConstPtrs(mf), 
                                varnames, geom, time, dt, ref_ratio);
#endif

    Vector<int> level_steps(nlevs, 0);
    WriteMultiLevelPlotfile("plt00000", nlevs, amrex::GetVecOfConstPtrs(mf),
                            varnames, geom, time, level_steps, ref_ratio);
   
    Vector<std::string> particle_realnames;
    for (int i = 0; i < NStructReal + NArrayReal; ++i)
    {
        particle_realnames.push_back("particle_real_component_" + std::to_string(i));
    }

    Vector<std::string> particle_intnames;
    for (int i = 0; i < NStructInt + NArrayInt; ++i)
    {
        particle_intnames.push_back("particle_int_component_" + std::to_string(i));
    }
    
    myPC.Checkpoint("plt00000", "particle0", false, particle_realnames, particle_intnames);
    }
    
    amrex::Finalize();
}
