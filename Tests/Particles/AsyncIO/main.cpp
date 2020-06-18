#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

static constexpr int NSR = 1 + AMREX_SPACEDIM;
static constexpr int NSI = 0;
static constexpr int NAR = 0;
static constexpr int NAI = 0;

struct TestParams {
  int nx;
  int ny;
  int nz;
  int max_grid_size;
  int nppc;
  int nlevs;
  bool verbose;
};

class MyParticleContainer
    : public amrex::ParticleContainer<NSR, NSI, NAR, NAI>
{

public:

    MyParticleContainer (const Vector<amrex::Geometry>            & a_geom,
                         const Vector<amrex::DistributionMapping> & a_dmap,
                         const Vector<amrex::BoxArray>            & a_ba,
                         const Vector<amrex::IntVect>             & a_rr)
        : amrex::ParticleContainer<NSR, NSI, NAR, NAI>(a_geom, a_dmap, a_ba, a_rr)
    {}

    std::size_t PSizeInFile (const Vector<int>& wrc, const Vector<int>& wic) const
    {
        std::size_t rsize = sizeof(ParticleReal)*std::accumulate(wrc.begin(), wrc.end(), 0);
        std::size_t isize = sizeof(int)*std::accumulate(wic.begin(), wic.end(), 0);
        return rsize + isize + AMREX_SPACEDIM*sizeof(ParticleReal) + 2*sizeof(int);
    }

    void WritePlotFileAsync (const std::string& dir, const std::string& name) const
    {
        Vector<int> wrc = {1, 1, 1, 1};
        Vector<int> wic;

        Vector<std::string> rcn = {"weight", "vx", "vy", "vz"};
        Vector<std::string> icn;
        WriteBinaryParticleDataAsync(dir, name, wrc, wic, rcn, icn);
    }

    void WriteBinaryParticleDataAsync (const std::string& dir, const std::string& name,
                                       const Vector<int>& write_real_comp,
                                       const Vector<int>& write_int_comp,
                                       const Vector<std::string>& real_comp_names,
                                       const Vector<std::string>& int_comp_names) const
    {
        BL_PROFILE("writeBinaryParticleDataAsync");

        AMREX_ASSERT(OK());

        AMREX_ASSERT(sizeof(typename ParticleType::RealType) == 4 ||
                     sizeof(typename ParticleType::RealType) == 8);

        const int MyProc = ParallelDescriptor::MyProc();
        const int NProcs = ParallelDescriptor::NProcs();

        AMREX_ALWAYS_ASSERT(real_comp_names.size() == NumRealComps() + NStructReal);
        AMREX_ALWAYS_ASSERT( int_comp_names.size() == NumIntComps() + NStructInt);

        std::string pdir = dir;
        if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
        pdir += name;

        if ( ! levelDirectoriesCreated)
        {
            if (ParallelDescriptor::IOProcessor())
                if ( ! amrex::UtilCreateDirectory(pdir, 0755))
                    amrex::CreateDirectoryFailed(pdir);
            ParallelDescriptor::Barrier();
        }

        Vector<Vector<Long> > np_per_grid(finestLevel()+1);;
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            np_per_grid[lev] = NumberOfParticlesInGrid(lev);
        }

        Vector<int64_t> np_on_rank(NProcs,-1L);
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
            {
                int rank = ParticleDistributionMap(lev)[mfi.index()];
                np_on_rank[rank] += np_per_grid[lev][mfi.index()];
            }
        }

        std::size_t psize = PSizeInFile(write_real_comp, write_int_comp);
        Vector<int64_t> rank_start_offset(NProcs);
        for (int ip = 0; ip < NProcs; ++ip) {
            auto info = AsyncOut::GetWriteInfo(ip);
            if (info.ispot == 0) {
                rank_start_offset[ip] = 0;
            } else {
                rank_start_offset[ip] = rank_start_offset[ip-1] + np_on_rank[ip-1]*psize;
            }
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::string HdrFileName = pdir;
            std::ofstream HdrFile;

            if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
                HdrFileName += '/';

            HdrFileName += "Header";
            HdrFileNamePrePost = HdrFileName;

            HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

            if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

            for (int lev = 0; lev <= finestLevel(); lev++)
            {
                Vector<int64_t> grid_offset(NProcs, 0);
                for (int k = 0; k < ParticleBoxArray(lev).size(); ++k)
                {
                    int rank = ParticleDistributionMap(lev)[k];
                    auto info = AsyncOut::GetWriteInfo(rank);
                    auto which = info.ifile;
                    auto count = np_per_grid[lev][k];
                    auto offset = grid_offset[rank] + rank_start_offset[rank];
                    HdrFile << which << ' ' << count << ' ' << offset << '\n';
                    grid_offset[rank] += count*psize;
                }
            }

            HdrFile.flush();
            HdrFile.close();
            if ( ! HdrFile.good())
            {
                amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
            }
        }

        using PinnedPTile = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
                                         PinnedArenaAllocator>;
        auto myptiles = std::make_shared<Vector<PinnedPTile> >();
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            for (MFIter mfi = MakeMFIter(lev, false); mfi.isValid(); ++mfi)
            {
                myptiles->emplace_back();
                auto& new_ptile = myptiles->back();

                if (np_per_grid[lev][mfi.index()] > 0)
                {
                    const auto& ptile = ParticlesAt(lev, mfi);
                    new_ptile.resize(ptile.numParticles());
                    amrex::copyParticles(new_ptile, ptile);
                }
            }
        }
    }
};

void test_async_io(TestParams& parms)
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

    IntVect domain_lo(D_DECL(0 , 0, 0));
    IntVect domain_hi(D_DECL(parms.nx - 1, parms.ny - 1, parms.nz-1));
    const Box domain(domain_lo, domain_hi);

    // Define the refinement ratio
    Vector<IntVect> rr(nlevs-1);
    for (int lev = 1; lev < nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));

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
        int n_fine = parms.nx*rr[0][0];
        IntVect refined_lo(D_DECL(n_fine/4,n_fine/4,n_fine/4));
        IntVect refined_hi(D_DECL(3*n_fine/4-1,3*n_fine/4-1,3*n_fine/4-1));

        // Build a box for the level 1 domain
        Box refined_patch(refined_lo, refined_hi);
        ba[1].define(refined_patch);
    }

    // break the BoxArrays at both levels into max_grid_size^3 boxes
    for (int lev = 0; lev < nlevs; lev++) {
        ba[lev].maxSize(parms.max_grid_size);
    }

    Vector<DistributionMapping> dmap(nlevs);

    Vector<std::unique_ptr<MultiFab> > partMF(nlevs);
    Vector<std::unique_ptr<MultiFab> > density(nlevs);
    Vector<std::unique_ptr<MultiFab> > acceleration(nlevs);
    for (int lev = 0; lev < nlevs; lev++) {
        dmap[lev] = DistributionMapping{ba[lev]};
        density[lev].reset(new MultiFab(ba[lev], dmap[lev], 1, 0));
        density[lev]->setVal(0.0);
        acceleration[lev].reset(new MultiFab(ba[lev], dmap[lev], 3, 1));
        acceleration[lev]->setVal(5.0, 1);
    }

    MyParticleContainer myPC(geom, dmap, ba, rr);
    myPC.SetVerbose(false);

    int num_particles = parms.nppc * AMREX_D_TERM(parms.nx, * parms.ny, * parms.nz);
    bool serialize = true;
    int iseed = 451;
    Real mass = 10.0;
    MyParticleContainer::ParticleInitData pdata = {mass};

    myPC.InitRandom(num_particles, iseed, pdata, serialize);

    myPC.AssignDensity(0, partMF, 0, 1, nlevs-1);

    for (int lev = 0; lev < nlevs; ++lev) {
        MultiFab::Copy(*density[lev], *partMF[lev], 0, 0, 1, 0);
    }

    Vector<std::string> varnames;
    varnames.push_back("density");

    Vector<std::string> particle_varnames;
    particle_varnames.push_back("mass");

    Vector<int> level_steps;
    level_steps.push_back(0);
    level_steps.push_back(0);

    int output_levs = nlevs;

    Vector<const MultiFab*> outputMF(output_levs);
    Vector<IntVect> outputRR(output_levs);
    for (int lev = 0; lev < output_levs; ++lev) {
        outputMF[lev] = density[lev].get();
        outputRR[lev] = IntVect(D_DECL(2, 2, 2));
    }

    WriteMultiLevelPlotfile("plt00000", output_levs, outputMF,
                            varnames, geom, 0.0, level_steps, outputRR);
    myPC.Checkpoint("plt00000", "particle0");
    myPC.WritePlotFileAsync("plt00000", "particle1");
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

  test_async_io(parms);

  amrex::Finalize();
}
