#include <iostream>

#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_AmrParticles.H>
#include <AMReX_PlotFileUtil.H>

#include <thread>

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
        const int IOProcNumber = NProcs - 1;

        AMREX_ALWAYS_ASSERT(real_comp_names.size() == NumRealComps() + NStructReal);
        AMREX_ALWAYS_ASSERT( int_comp_names.size() == NumIntComps() + NStructInt);

        std::string pdir = dir;
        if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
        pdir += name;

        if ( ! levelDirectoriesCreated)
        {
            if (MyProc == IOProcNumber)
                if ( ! amrex::UtilCreateDirectory(pdir, 0755))
                    amrex::CreateDirectoryFailed(pdir);
            ParallelDescriptor::Barrier();
        }

        Long total_np = 0;
        Vector<Long> np_per_level(finestLevel()+1);
        Vector<Vector<Long> > np_per_grid(finestLevel()+1);
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            np_per_grid[lev] = NumberOfParticlesInGrid(lev);
            np_per_level[lev] = std::accumulate(np_per_grid[lev].begin(),
                                                np_per_grid[lev].end(), 0);
            total_np += np_per_level[lev];
        }


        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            std::string LevelDir = pdir;
            bool gotsome = np_per_level[lev];

            if (gotsome)
            {
                if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';

                LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);

                if ( ! levelDirectoriesCreated) {
                    if (MyProc == IOProcNumber)
                        if ( ! amrex::UtilCreateDirectory(LevelDir, 0755))
                            amrex::CreateDirectoryFailed(LevelDir);
                    ParallelDescriptor::Barrier();
                }
            }

            // Write out the header for each particle
            if (gotsome and (MyProc == IOProcNumber))
            {
                std::string HeaderFileName = LevelDir;
                HeaderFileName += "/Particle_H";
                std::ofstream ParticleHeader(HeaderFileName);

                ParticleBoxArray(lev).writeOn(ParticleHeader);
                ParticleHeader << '\n';

                ParticleHeader.flush();
                ParticleHeader.close();
            }
        }

        int maxnextid = ParticleType::NextID();
        ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);

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

        // make tmp particle tiles in pinned memory to write
        using PinnedPTile = ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
                                         PinnedArenaAllocator>;
        auto myptiles = std::make_shared<Vector<std::map<std::pair<int, int>,PinnedPTile> > >();
        myptiles->resize(finestLevel()+1);
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                auto& new_ptile = (*myptiles)[lev][std::make_pair(mfi.index(),
                                                                  mfi.LocalTileIndex())];

                if (np_per_grid[lev][mfi.index()] > 0)
                {
                    const auto& ptile = ParticlesAt(lev, mfi);
                    new_ptile.resize(ptile.numParticles());
                    amrex::copyParticles(new_ptile, ptile);
                }
            }
        }

        auto wrc = std::make_shared<Vector<int> >(write_real_comp);
        auto wic = std::make_shared<Vector<int> >(write_int_comp);
        auto rcn = std::make_shared<Vector<std::string> >(real_comp_names);
        auto icn = std::make_shared<Vector<std::string> >(int_comp_names);

        int finest_level = finestLevel();
        Vector<BoxArray> bas;
        Vector<DistributionMapping> dms;
        for (int lev = 0; lev <= finestLevel(); lev++)
        {
            bas.push_back(ParticleBoxArray(lev));
            dms.push_back(ParticleDistributionMap(lev));
        }

        int nrc = NumRealComps();
        int nic = NumIntComps();

        auto RD = ParticleRealDescriptor;

        AsyncOut::Submit([=] ()
        {
            if (MyProc == IOProcNumber)
            {
                std::string HdrFileName = pdir;
                std::ofstream HdrFile;

                if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
                    HdrFileName += '/';

                HdrFileName += "Header";

                HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

                if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

                if (sizeof(typename ParticleType::RealType) == 4)
                {
                    HdrFile << ParticleType::Version() << "_single" << '\n';
                }
                else
                {
                    HdrFile << ParticleType::Version() << "_double" << '\n';
                }

                int num_output_real = 0;
                for (int i = 0; i < nrc + NStructReal; ++i)
                    if ((*wrc)[i]) ++num_output_real;

                int num_output_int = 0;
                for (int i = 0; i < nic + NStructInt; ++i)
                    if ((*wic)[i]) ++num_output_int;

                // AMREX_SPACEDIM and N for sanity checking.
                HdrFile << AMREX_SPACEDIM << '\n';

                // The number of extra real parameters
                HdrFile << num_output_real << '\n';

                // Real component names
                for (int i = 0; i < NStructReal + nrc; ++i )
                    if ((*wrc)[i]) HdrFile << (*rcn)[i] << '\n';

                // The number of extra int parameters
                HdrFile << num_output_int << '\n';

                // int component names
                for (int i = 0; i < NStructInt + nic; ++i )
                    if ((*wic)[i]) HdrFile << (*icn)[i] << '\n';

                bool is_checkpoint = true; // legacy
                HdrFile << is_checkpoint << '\n';

                // The total number of particles.
                HdrFile << total_np << '\n';

                // The value of nextid that we need to restore on restart.
                HdrFile << maxnextid << '\n';

                // Then the finest level of the AMR hierarchy.
                HdrFile << finest_level << '\n';

                // Then the number of grids at each level.
                for (int lev = 0; lev <= finest_level; lev++)
                    HdrFile << dms[lev].size() << '\n';

                for (int lev = 0; lev <= finest_level; lev++)
                {
                    Vector<int64_t> grid_offset(NProcs, 0);
                    for (int k = 0; k < bas[lev].size(); ++k)
                    {
                        int rank = dms[lev][k];
                        auto info = AsyncOut::GetWriteInfo(rank);
                        HdrFile << info.ifile << ' '
                                << np_per_grid[lev][k] << ' '
                                << grid_offset[rank] + rank_start_offset[rank] << '\n';
                        grid_offset[rank] += np_per_grid[lev][k]*psize;
                    }
                }

                HdrFile.flush();
                HdrFile.close();
                if ( ! HdrFile.good())
                {
                    amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
                }
            }

            AsyncOut::Wait();  // Wait for my turn

            for (int lev = 0; lev <= finest_level; lev++)
            {
                // For a each grid, the tiles it contains
                std::map<int, Vector<int> > tile_map;

                for (const auto& kv : (*myptiles)[lev])
                {
                    const int grid = kv.first.first;
                    const int tile = kv.first.second;
                    tile_map[grid].push_back(tile);
                }

                std::string LevelDir = pdir;
                if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';
                LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);
                std::string filePrefix(LevelDir);
                filePrefix += '/';
                filePrefix += ParticleType::DataPrefix();
                auto info = AsyncOut::GetWriteInfo(MyProc);
                std::string file_name = amrex::Concatenate(filePrefix, info.ifile, 5);
                std::ofstream ofs;
                ofs.open(file_name.c_str(), (info.ispot == 0) ? (std::ios::binary | std::ios::trunc)
                         : (std::ios::binary | std::ios::app));

                for (int k = 0; k < bas[lev].size(); ++k)
                {
                    int rank = dms[lev][k];
                    if (rank != MyProc) continue;
                    const int grid = k;
                    if (np_per_grid[lev][grid] == 0) continue;

                    // First write out the integer data in binary.
                    int num_output_int = 0;
                    for (int i = 0; i < nic + NStructInt; ++i)
                        if ((*wic)[i]) ++num_output_int;

                    const int iChunkSize = 2 + num_output_int;
                    Vector<int> istuff(np_per_grid[lev][grid]*iChunkSize);
                    int* iptr = istuff.dataPtr();

                    for (unsigned i = 0; i < tile_map[grid].size(); i++) {
                        auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
                        const auto& pbox = (*myptiles)[lev][ptile_index];
                        for (int pindex = 0;
                             pindex < pbox.GetArrayOfStructs().numParticles(); ++pindex)
                        {
                            const auto& aos = pbox.GetArrayOfStructs();
                            const auto& p = aos[pindex];

                            // always write these
                            *iptr = p.id(); ++iptr;
                            *iptr = p.cpu(); ++iptr;

                            // optionally write these
                            for (int j = 0; j < NStructInt; j++)
                            {
                                if (write_int_comp[j])
                                {
                                    *iptr = p.idata(j);
                                    ++iptr;
                                }
                            }

                            const auto& soa  = pbox.GetStructOfArrays();
                            for (int j = 0; j < nic; j++)
                            {
                                if (write_int_comp[NStructInt+j])
                                {
                                    *iptr = soa.GetIntData(j)[pindex];
                                    ++iptr;
                                }
                            }
                        }
                    }

                    writeIntData(istuff.dataPtr(), istuff.size(), ofs);
                    ofs.flush();  // Some systems require this flush() (probably due to a bug)

                    // Write the Real data in binary.
                    int num_output_real = 0;
                    for (int i = 0; i < nrc + NStructReal; ++i)
                        if ((*wrc)[i]) ++num_output_real;

                    const int rChunkSize = AMREX_SPACEDIM + num_output_real;
                    Vector<typename ParticleType::RealType> rstuff(np_per_grid[lev][grid]*rChunkSize);
                    typename ParticleType::RealType* rptr = rstuff.dataPtr();

                    for (unsigned i = 0; i < tile_map[grid].size(); i++) {
                        auto ptile_index = std::make_pair(grid, tile_map[grid][i]);
                        const auto& pbox = (*myptiles)[lev][ptile_index];
                        for (int pindex = 0;
                             pindex < pbox.GetArrayOfStructs().numParticles(); ++pindex)
                        {
                            const auto& aos = pbox.GetArrayOfStructs();
                            const auto& p = aos[pindex];

                            // always write these
                            for (int j = 0; j < AMREX_SPACEDIM; j++) rptr[j] = p.pos(j);
                            rptr += AMREX_SPACEDIM;

                            // optionally write these
                            for (int j = 0; j < NStructReal; j++)
                            {
                                if (write_real_comp[j])
                                {
                                    *rptr = p.rdata(j);
                                    ++rptr;
                                }
                            }

                            const auto& soa  = pbox.GetStructOfArrays();
                            for (int j = 0; j < nrc; j++)
                            {
                                if (write_real_comp[NStructReal+j])
                                {
                                    *rptr = (typename ParticleType::RealType) soa.GetRealData(j)[pindex];
                                    ++rptr;
                                }
                            }
                        }
                    }

                    if (sizeof(typename ParticleType::RealType) == 4) {
                        writeFloatData((float*) rstuff.dataPtr(), rstuff.size(), ofs, RD);
                    }
                    else if (sizeof(typename ParticleType::RealType) == 8) {
                        writeDoubleData((double*) rstuff.dataPtr(), rstuff.size(), ofs, RD);
                    }

                    ofs.flush();  // Some systems require this flush() (probably due to a bug)
                }
            }
            AsyncOut::Notify();  // Notify others I am done
        });
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
