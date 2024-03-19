#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

#if !defined(AMREX_PARTICLES) || !defined(AMREX_USE_SENSEI_INSITU)
#error Incompatible AMReX library configuration! This tutorial requires AMREX_PARTICLES and AMREX_USE_SENSEI_INSITU
#endif
#include <AMReX_ParticleInSituBridge.H>

using namespace amrex;

static constexpr int NR = 7;
static constexpr int NI = 4;

int num_runtime_real = 0;
int num_runtime_int = 0;

bool remove_negative = true;

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
    : public amrex::ParticleContainerPureSoA<NR, NI>
{

public:

    TestParticleContainer (const Vector<amrex::Geometry>            & a_geom,
                           const Vector<amrex::DistributionMapping> & a_dmap,
                           const Vector<amrex::BoxArray>            & a_ba,
                           const Vector<amrex::IntVect>             & a_rr)
        : amrex::ParticleContainerPureSoA<NR, NI>(a_geom, a_dmap, a_ba, a_rr)
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

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell)
    {
        BL_PROFILE("InitParticles");

        const int lev = 0;  // only add particles on level 0
        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();

        const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                         *a_num_particles_per_cell[1],
                                         *a_num_particles_per_cell[2]);

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();

            std::array<Gpu::HostVector<ParticleReal>, NR> host_real;
            std::array<Gpu::HostVector<int>, NI> host_int;

            std::vector<Gpu::HostVector<ParticleReal> > host_runtime_real(NumRuntimeRealComps());
            std::vector<Gpu::HostVector<int> > host_runtime_int(NumRuntimeIntComps());

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    Real r[3];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                    amrex::Long id = ParticleType::NextID();

                    host_int[0].push_back(static_cast<int>(id));
                    host_int[1].push_back(ParallelDescriptor::MyProc());
                    host_real[0].push_back(static_cast<ParticleReal> (plo[0] + (iv[0] + r[0])*dx[0]));
#if AMREX_SPACEDIM > 1
                    host_real[1].push_back(static_cast<ParticleReal> (plo[1] + (iv[1] + r[1])*dx[1]));
#endif
#if AMREX_SPACEDIM > 2
                    host_real[2].push_back(static_cast<ParticleReal> (plo[2] + (iv[2] + r[2])*dx[2]));
#endif

                    for (int i = AMREX_SPACEDIM; i < NR; ++i)
                        host_real[i].push_back(static_cast<ParticleReal>(id));
                    for (int i = 2; i < NI; ++i)
                        host_int[i].push_back(static_cast<int>(id));
                    for (int i = 0; i < NumRuntimeRealComps(); ++i)
                        host_runtime_real[i].push_back(static_cast<ParticleReal>(id));
                    for (int i = 0; i < NumRuntimeIntComps(); ++i)
                        host_runtime_int[i].push_back(static_cast<int>(id));
                }
            }

            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            auto old_size = particle_tile.size();
            auto new_size = old_size + host_real[0].size();
            particle_tile.resize(new_size);

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NR; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_real[i].begin(),
                               host_real[i].end(),
                               soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NI; ++i)
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
                               soa.GetRealData(NR+i).begin() + old_size);
            }

            for (int i = 0; i < NumRuntimeIntComps(); ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice,
                               host_runtime_int[i].begin(),
                               host_runtime_int[i].end(),
                               soa.GetIntData(NI+i).begin() + old_size);
            }

            Gpu::streamSynchronize();
        }

        Redistribute();
    }
};

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int num_ppc;
    int is_periodic;
    int nlevs;
};

void testRedistribute();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running redistribute test \n";
    testRedistribute();

    amrex::Finalize();
}

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
    pp.get("is_periodic", params.is_periodic);
    pp.get("nlevs", params.nlevs);
    pp.query("num_runtime_real", num_runtime_real);
    pp.query("num_runtime_int", num_runtime_int);
}

void testRedistribute ()
{
    BL_PROFILE("testSENSEI_insitu");
    TestParams params;
    get_test_params(params, "insitu");

    int is_per[BL_SPACEDIM];
    for (int & d : is_per)
        d = params.is_periodic;

    Vector<IntVect> rr(params.nlevs-1);
    for (int lev = 1; lev < params.nlevs; lev++)
        rr[lev-1] = IntVect(AMREX_D_DECL(2,2,2));

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box base_domain(domain_lo, domain_hi);

    Vector<Geometry> geom(params.nlevs);
    geom[0].define(base_domain, &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < params.nlevs; lev++) {
        geom[lev].define(amrex::refine(geom[lev-1].Domain(), rr[lev-1]),
                         &real_box, CoordSys::cartesian, is_per);
    }

    Vector<BoxArray> ba(params.nlevs);
    Vector<DistributionMapping> dm(params.nlevs);
    auto lo = IntVect(AMREX_D_DECL(0, 0, 0));
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

    int npc = params.num_ppc;
    auto nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    amrex::Print() << "About to initialize particles \n";

    pc.InitParticles(nppc);

    auto *insitu_bridge = new ParticleInSituBridge;

    if (insitu_bridge->initialize()) {
        amrex::ErrorStream() << "Failed to initialize the in situ bridge." << '\n';
        amrex::Abort();
    }

    // define specifications for fields on particles
    std::map<std::string, std::vector<int>> rStructs;
    std::map<std::string, int> iStructs;
    std::map<std::string, std::vector<int>> rArrays;
    std::map<std::string, int> iArrays;

    if (insitu_bridge->update(0.0, 0, &pc, rStructs)) {
        amrex::ErrorStream() << "Failed to update the in situ bridge." << '\n';
        amrex::Abort();
    }

    if (insitu_bridge->finalize()) {
        amrex::ErrorStream() << "Failed to finalize the in situ bridge." << '\n';
    }

    delete insitu_bridge;

    // the way this test is set up, if we make it here we pass
    amrex::Print() << "pass \n";
}
