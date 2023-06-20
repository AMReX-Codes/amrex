#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_Particles.H>

#include <unistd.h>
#include <cstdio>

using namespace amrex;

void set_grids_nested (Vector<Box>& domains,
                       Vector<BoxArray>& grids,
                       Vector<IntVect>& ref_ratio);
void set_grids_file (Vector<Box>& domains,
                     Vector<BoxArray>& grids,
                     Vector<IntVect>& ref_ratio);
void test ();
BoxArray readBoxList (const std::string& file, Box& domain);

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    test();
    amrex::Finalize();
}

void test ()
{
    const int nghost = 0;
    int ncells, max_grid_size, ncomp, nlevs, nppc;
    int restart_check = 0, nplotfile = 1, nparticlefile = 1, sleeptime = 0;
    int grids_from_file = 0;
    std::string compression = "None@0";
    std::string directory = "";

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("ncomp", ncomp);
    pp.get("nlevs", nlevs);
    pp.get("nppc", nppc);
    pp.query("hdf5compression", compression);
    pp.query("nplotfile", nplotfile);
    pp.query("nparticlefile", nparticlefile);
    pp.query("sleeptime", sleeptime);
    pp.query("restart_check", restart_check);
    pp.query("grids_from_file", grids_from_file);
    pp.query("directory", directory);

    if (directory != "" && directory.back() != '/') {
        // Include separator if one was not provided
        directory += "/";
    }

    Vector<Box> domains;
    Vector<BoxArray> ba;
    Vector<IntVect> ref_ratio;

    if (grids_from_file)
    {
        set_grids_file(domains, ba, ref_ratio);
    }
    else
    {
        set_grids_nested(domains, ba, ref_ratio);
    }

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }

    // This sets the boundary conditions to be doubly or triply periodic
    int is_per[AMREX_SPACEDIM];
    for (int i = 0; i < AMREX_SPACEDIM; i++)
        is_per[i] = 1;

    // This defines a Geometry object for each level
    Vector<Geometry> geom(nlevs);
    geom[0].define(domains[0], &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; lev++) {
        geom[lev].define(domains[lev], &real_box, CoordSys::cartesian, is_per);
    }

    Vector<DistributionMapping> dmap(nlevs);

    Vector<std::unique_ptr<MultiFab> > mf(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        mf[lev] = std::make_unique<MultiFab>(ba[lev], dmap[lev], ncomp, nghost);
        mf[lev]->setVal(lev);
    }

    // these don't really matter, make something up
    const Real time = 0.0;
    /* const Real dt = 0.0; */

    Vector<std::string> varnames;
    for (int i = 0; i < ncomp; ++i)
    {
        varnames.push_back("component_" + std::to_string(i));
    }

    Vector<int> level_steps(nlevs, 0);

    /* if (compression.compare("None@0") != 0) */
    /*     std::cout << "Compression: " << compression << std::endl; */

    char fname[512];
    for (int ts = 0; ts < nplotfile; ts++) {
        std::snprintf(fname, sizeof fname, "%splt%05d", directory.c_str(), ts);

        // Fake computation
        if (ts > 0 && sleeptime > 0) {
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "Sleep for " << sleeptime << " seconds." << std::endl;
                fflush(stdout);
            }
            sleep(sleeptime);
        }

        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Writing plot file [" << fname << ".h5] ...";
            fflush(stdout);
        }
#ifdef AMREX_USE_HDF5
#if (defined AMREX_USE_HDF5_ZFP) || (defined AMREX_USE_HDF5_SZ)
        WriteMultiLevelPlotfileHDF5SingleDset(fname, nlevs, amrex::GetVecOfConstPtrs(mf), varnames,
                                             geom, time, level_steps, ref_ratio, compression);
        /* WriteMultiLevelPlotfileHDF5MultiDset(fname, nlevs, amrex::GetVecOfConstPtrs(mf), varnames, */
        /*                                      geom, time, level_steps, ref_ratio, compression); */
#else
        WriteMultiLevelPlotfileHDF5SingleDset(fname, nlevs, amrex::GetVecOfConstPtrs(mf), varnames,
                                              geom, time, level_steps, ref_ratio);
#endif
#else
        WriteMultiLevelPlotfile(fname, nlevs, amrex::GetVecOfConstPtrs(mf),
                                varnames, geom, time, level_steps, ref_ratio);
#endif
        if (ParallelDescriptor::IOProcessor())
            std::cout << " done" << std::endl;
    }

    /* ParallelDescriptor::Barrier(); */

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

    if (nparticlefile > 0) {
        if (ParallelDescriptor::IOProcessor()) {
            std::cout << "Init particles ...";
            fflush(stdout);
        }

        myPC.InitRandom(num_particles, iseed, pdata, serialize);

        if (ParallelDescriptor::IOProcessor())
            std::cout << " done" << std::endl;

        Vector<std::string> particle_realnames;
        for (int i = 0; i < NStructReal + NArrayReal; ++i)
            particle_realnames.push_back("particle_real_component_" + std::to_string(i));

        Vector<std::string> particle_intnames;
        for (int i = 0; i < NStructInt + NArrayInt; ++i)
            particle_intnames.push_back("particle_int_component_" + std::to_string(i));

        for (int ts = 0; ts < nparticlefile; ts++) {
            std::snprintf(fname, sizeof fname, "%splt%05d", directory.c_str(), ts);

            // Fake computation
            if (ts > 0 && sleeptime > 0) {
                if (ParallelDescriptor::IOProcessor()) {
                    std::cout << "Sleep for " << sleeptime << " seconds." << std::endl;
                    fflush(stdout);
                }
                sleep(sleeptime);
            }

            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "Writing particle file [" << fname << "/particle0.h5] ...";
                fflush(stdout);
            }

#ifdef AMREX_USE_HDF5
            myPC.CheckpointHDF5(fname, "particle0", false, particle_realnames, particle_intnames, compression);
#else
            myPC.Checkpoint(fname, "particle0", false, particle_realnames, particle_intnames);
            /* myPC.WriteAsciiFile("particle0_ascii"); */
#endif
            if (ParallelDescriptor::IOProcessor())
                std::cout << " done" << std::endl;
        }
    }

    ParallelDescriptor::Barrier();

    char directory_path[512];
    if (restart_check && nparticlefile > 0)
    {
        MyPC newPC(geom, dmap, ba, ref_ratio);
#ifdef AMREX_USE_HDF5
        std::snprintf(directory_path, sizeof directory_path, "%s%s", directory.c_str(), "plt00000/particle0");
        newPC.RestartHDF5(directory_path, "particle0");
#else
        std::snprintf(directory_path, sizeof directory_path, "%s%s", directory.c_str(), "plt00000");
        newPC.Restart(directory_path, "particle0");
#endif

        using PType = typename MyPC::SuperParticleType;

        for (int icomp=0; icomp<NStructReal+NArrayReal+NStructInt+NArrayInt; ++icomp)
        {
            auto sm_new = amrex::ReduceSum(newPC,
                [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                {
                    return p.rdata(1);
                });

            auto sm_old = amrex::ReduceSum(myPC,
                [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                {
                    return p.rdata(1);
                });

            ParallelDescriptor::ReduceRealSum(sm_new);
            ParallelDescriptor::ReduceRealSum(sm_old);

            AMREX_ALWAYS_ASSERT(sm_old == sm_new);
        }
    }
}

BoxArray
readBoxList (const std::string& file, Box& domain)
{
    BoxArray retval;

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(file, fileCharPtr);
    std::istringstream is(fileCharPtr.data());

    is >> domain;
    is.ignore(1000,'\n');
    retval.readFrom(is);

    return retval;
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
    for (int lev = 1; lev < nlevs; lev++)
        ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));

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

void set_grids_file (Vector<Box>& domains,
                     Vector<BoxArray>& grids,
                     Vector<IntVect>& ref_ratio)
{
    int ncells, max_grid_size, nlevs;
    Vector<int> ref_ratio_file;

    ParmParse pp;
    pp.get("ncells", ncells);
    pp.get("max_grid_size", max_grid_size);
    pp.get("nlevs", nlevs);

    ref_ratio_file.resize(nlevs-1);
    pp.getarr("ref_ratio_file", ref_ratio_file);

    domains.resize(nlevs);
    grids.resize(nlevs);
    ref_ratio.resize(nlevs-1);

    for (int lev=0; lev < nlevs-1; ++lev)
    {
        ref_ratio[lev] = IntVect(AMREX_D_DECL(ref_ratio_file[lev],
                                              ref_ratio_file[lev],
                                              ref_ratio_file[lev]));
    }

    for (int lev=0; lev < nlevs; ++lev)
    {
        Box domain;
        grids[lev] = readBoxList("grids/Level_"+std::to_string(lev), domain);
        domains[lev] = domain;
    }
}
