#include <AMReX.H>
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>

#include <algorithm>
#include <cstdint>
#include <numeric>
#include <random>

using namespace amrex;

static constexpr int NSR = 6;
static constexpr int NSI = 1;
static constexpr int NAR = 1;
static constexpr int NAI = 1;

int num_runtime_real = 0;
int num_runtime_int = 0;

void get_position_unit_cell (Real* r, const IntVect& nppc, int i_part)
{
    int nx = nppc[0];
#if AMREX_SPACEDIM > 1
    int ny = nppc[1];
#else
    int ny = 1;
#endif
#if AMREX_SPACEDIM > 2
    int nz = nppc[2];
#else
    int nz = 1;
#endif

    int ix_part = i_part/(ny * nz);
    int iy_part = (i_part % (ny * nz)) % ny;
    int iz_part = (i_part % (ny * nz)) / ny;

    r[0] = (0.5+ix_part)/nx;
    r[1] = (0.5+iy_part)/ny;
    r[2] = (0.5+iz_part)/nz;
}

class TestParticleContainer
    : public amrex::ParticleContainer<NSR, NSI, NAR, NAI>
{
public:

    TestParticleContainer (const Vector<amrex::Geometry>& a_geom,
                           const Vector<amrex::DistributionMapping>& a_dmap,
                           const Vector<amrex::BoxArray>& a_ba,
                           const Vector<amrex::IntVect>& a_rr)
        : amrex::ParticleContainer<NSR, NSI, NAR, NAI>(a_geom, a_dmap, a_ba, a_rr)
    {
        for (int i = 0; i < num_runtime_real; ++i)
        {
            AddRealComp(true);
        }
        for (int i = 0; i < num_runtime_int; ++i)
        {
            AddIntComp(true);
        }
    }

    void RedistributeGlobal ()
    {
        const int lev_min = 0;
        const int lev_max = finestLevel();
        const int nGrow = 0;
        const int local = 0;
        Redistribute(lev_min, lev_max, nGrow, local);
    }

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell)
    {
        BL_PROFILE("InitParticles");

        const int lev = 0;
        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        const int num_ppc = AMREX_D_TERM(a_num_particles_per_cell[0],
                                         * a_num_particles_per_cell[1],
                                         * a_num_particles_per_cell[2]);

        for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box = mfi.tilebox();

            Gpu::HostVector<ParticleType> host_particles;
            std::array<Gpu::HostVector<ParticleReal>, NAR> host_real;
            std::array<Gpu::HostVector<int>, NAI> host_int;

            std::vector<Gpu::HostVector<ParticleReal> > host_runtime_real(NumRuntimeRealComps());
            std::vector<Gpu::HostVector<int> > host_runtime_int(NumRuntimeIntComps());

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part = 0; i_part < num_ppc; ++i_part)
                {
                    Real r[3];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                    ParticleType p;
                    p.id() = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    p.pos(0) = static_cast<ParticleReal>(plo[0] + (iv[0] + r[0])*dx[0]);
#if AMREX_SPACEDIM > 1
                    p.pos(1) = static_cast<ParticleReal>(plo[1] + (iv[1] + r[1])*dx[1]);
#endif
#if AMREX_SPACEDIM > 2
                    p.pos(2) = static_cast<ParticleReal>(plo[2] + (iv[2] + r[2])*dx[2]);
#endif

                    for (int i = 0; i < NSR; ++i) { p.rdata(i) = ParticleReal(p.id()); }
                    for (int i = 0; i < NSI; ++i) { p.idata(i) = int(p.id()); }

                    host_particles.push_back(p);
                    for (int i = 0; i < NAR; ++i) {
                        host_real[i].push_back(ParticleReal(p.id()));
                    }
                    for (int i = 0; i < NAI; ++i) {
                        host_int[i].push_back(int(p.id()));
                    }
                    for (int i = 0; i < NumRuntimeRealComps(); ++i) {
                        host_runtime_real[i].push_back(ParticleReal(p.id()));
                    }
                    for (int i = 0; i < NumRuntimeIntComps(); ++i) {
                        host_runtime_int[i].push_back(int(p.id()));
                    }
                }
            }

            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);

            Gpu::copyAsync(Gpu::hostToDevice,
                           host_particles.begin(),
                           host_particles.end(),
                           particle_tile.GetArrayOfStructs().begin() + old_size);

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_real[i].begin(),
                               host_real[i].end(),
                               soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_int[i].begin(),
                               host_int[i].end(),
                               soa.GetIntData(i).begin() + old_size);
            }

            for (int i = 0; i < NumRuntimeRealComps(); ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_runtime_real[i].begin(),
                               host_runtime_real[i].end(),
                               soa.GetRealData(NAR+i).begin() + old_size);
            }

            for (int i = 0; i < NumRuntimeIntComps(); ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_runtime_int[i].begin(),
                               host_runtime_int[i].end(),
                               soa.GetIntData(NAI+i).begin() + old_size);
            }

            Gpu::streamSynchronize();
        }

        RedistributeGlobal();
    }

    void checkAnswer () const
    {
        BL_PROFILE("TestParticleContainer::checkAnswer");

        AMREX_ALWAYS_ASSERT(OK());

        int num_rr = NumRuntimeRealComps();
        int num_ii = NumRuntimeIntComps();

        for (int lev = 0; lev <= finestLevel(); ++lev)
        {
            const auto& plev = GetParticles(lev);
            for (MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                int gid = mfi.index();
                int tid = mfi.LocalTileIndex();
                const auto& ptile = plev.at(std::make_pair(gid, tid));
                const auto& ptd = ptile.getConstParticleTileData();
                const size_t np = ptile.numParticles();

                AMREX_FOR_1D(np, i,
                {
                    for (int j = 0; j < NSR; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_aos[i].rdata(j) == ptd.m_aos[i].id());
                    }
                    for (int j = 0; j < NSI; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_aos[i].idata(j) == ptd.m_aos[i].id());
                    }
                    if constexpr (NAR > 0) {
                        for (int j = 0; j < NAR; ++j)
                        {
                            AMREX_ALWAYS_ASSERT(ptd.m_rdata[j][i] == ptd.m_aos[i].id());
                        }
                    }
                    if constexpr (NAI > 0) {
                        for (int j = 0; j < NAI; ++j)
                        {
                            AMREX_ALWAYS_ASSERT(ptd.m_idata[j][i] == ptd.m_aos[i].id());
                        }
                    }
                    for (int j = 0; j < num_rr; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_runtime_rdata[j][i] == ptd.m_aos[i].id());
                    }
                    for (int j = 0; j < num_ii; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_runtime_idata[j][i] == ptd.m_aos[i].id());
                    }
                });
            }
        }
    }
};

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int num_ppc;
    int is_periodic;
    int nsteps;
    int nlevs;
    int sort;
    int stable_redistribute = 0;
    int random_seed = 8675309;
    int check_answer_each_step = 1;
};

auto makeRandomPMap (int nboxes, int nprocs, std::uint32_t seed) -> Vector<int>
{
    Vector<int> pmap(nboxes);
    for (int i = 0; i < nboxes; ++i) {
        pmap[i] = i % nprocs;
    }

    std::mt19937 gen(seed);
    std::shuffle(pmap.begin(), pmap.end(), gen);

    return pmap;
}

void testRedistributeGlobalDM ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running global redistribute DistributionMap shuffle test\n";
    testRedistributeGlobalDM();

    amrex::Finalize();
}

void get_test_params (TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
    pp.get("is_periodic", params.is_periodic);
    pp.get("nsteps", params.nsteps);
    pp.get("nlevs", params.nlevs);
    pp.query("num_runtime_real", num_runtime_real);
    pp.query("num_runtime_int", num_runtime_int);
    pp.query("stable_redistribute", params.stable_redistribute);
    pp.query("random_seed", params.random_seed);
    pp.query("check_answer_each_step", params.check_answer_each_step);

    params.sort = 0;
    pp.query("sort", params.sort);
}

void testRedistributeGlobalDM ()
{
    BL_PROFILE("testRedistributeGlobalDM");
    TestParams params;
    get_test_params(params, "redistribute_global_dm");

    int is_per[] = {AMREX_D_DECL(params.is_periodic,
                                 params.is_periodic,
                                 params.is_periodic)};

    Vector<IntVect> rr(params.nlevs-1);
    for (int lev = 1; lev < params.nlevs; ++lev) {
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));
    }

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; ++n)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box base_domain(domain_lo, domain_hi);

    Vector<Geometry> geom(params.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < params.nlevs; ++lev) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(params.nlevs);
    Vector<DistributionMapping> dm(params.nlevs);
    IntVect lo(0);
    IntVect size = params.size;
    for (int lev = 0; lev < params.nlevs; ++lev)
    {
        ba[lev].define(Box(lo, lo+params.size-1));
        ba[lev].maxSize(params.max_grid_size);
        dm[lev].define(ba[lev]);
        lo += size/2;
        size *= 2;
    }

    TestParticleContainer pc(geom, dm, ba, rr);
    pc.setStableRedistribute(params.stable_redistribute);

    IntVect nppc(params.num_ppc);

    amrex::Print() << "About to initialize particles\n";

    pc.InitParticles(nppc);
    pc.checkAnswer();

    auto np_old = pc.TotalNumberOfParticles();
    const int nprocs = ParallelDescriptor::NProcs();

    amrex::Print() << "Benchmark setup: " << ba[0].size() << " boxes on level 0 across "
                   << nprocs << " MPI ranks\n";

    if (params.sort) { pc.SortParticlesByCell(); }

    Real total_dm_time = Real(0.0);
    Real total_redistribute_time = Real(0.0);

    for (int i = 0; i < params.nsteps; ++i)
    {
        const auto dm_start = amrex::second();
        for (int lev = 0; lev < params.nlevs; ++lev)
        {
            auto pmap = makeRandomPMap(ba[lev].size(), nprocs,
                                       static_cast<std::uint32_t>(params.random_seed + 7919*i + 101*lev));
            DistributionMapping new_dm;
            new_dm.define(pmap);
            pc.SetParticleDistributionMap(lev, new_dm);
        }
        total_dm_time += amrex::second() - dm_start;

        ParallelDescriptor::Barrier();
        const auto redistribute_start = amrex::second();
        pc.RedistributeGlobal();
        total_redistribute_time += amrex::second() - redistribute_start;

        if (params.sort) { pc.SortParticlesByCell(); }
        if (params.check_answer_each_step) { pc.checkAnswer(); }
    }

    if (!params.check_answer_each_step) {
        pc.checkAnswer();
    }

    if (geom[0].isAllPeriodic()) {
        AMREX_ALWAYS_ASSERT(np_old == pc.TotalNumberOfParticles());
    }

    ParallelDescriptor::ReduceRealMax(total_dm_time, ParallelDescriptor::IOProcessorNumber());
    ParallelDescriptor::ReduceRealMax(total_redistribute_time, ParallelDescriptor::IOProcessorNumber());

    amrex::Print() << "Max DM shuffle time over all ranks: " << total_dm_time << " s\n";
    amrex::Print() << "Max redistribute time over all ranks: " << total_redistribute_time << " s\n";
    amrex::Print() << "Average redistribute time per step: "
                   << total_redistribute_time/static_cast<Real>(params.nsteps) << " s\n";
    amrex::Print() << "pass\n";
}
