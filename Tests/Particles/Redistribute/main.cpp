#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

using namespace amrex;

static constexpr int NSR = 4;
static constexpr int NSI = 3;
static constexpr int NAR = 2;
static constexpr int NAI = 1;

int num_runtime_real = 0;
int num_runtime_int = 0;

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
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

    TestParticleContainer (const Vector<amrex::Geometry>            & a_geom,
                           const Vector<amrex::DistributionMapping> & a_dmap,
                           const Vector<amrex::BoxArray>            & a_ba,
                           const Vector<amrex::IntVect>             & a_rr)
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

    void RedistributeLocal ()
    {
        const int lev_min = 0;
        const int lev_max = finestLevel();
        const int nGrow = 0;
        const int local = 1;
        Redistribute(lev_min, lev_max, nGrow, local);
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
        
        const int lev = 0;  // only add particles on level 0
        const Real* dx = Geom(lev).CellSize();
        const Real* plo = Geom(lev).ProbLo();
    
        const int num_ppc = AMREX_D_TERM( a_num_particles_per_cell[0],
                                         *a_num_particles_per_cell[1],
                                         *a_num_particles_per_cell[2]);

        for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
        {
            const Box& tile_box  = mfi.tilebox();

            Gpu::HostVector<ParticleType> host_particles;
            std::array<Gpu::HostVector<ParticleReal>, NAR> host_real;
            std::array<Gpu::HostVector<int>, NAI> host_int;

            std::vector<Gpu::HostVector<ParticleReal> > host_runtime_real(NumRuntimeRealComps());
            std::vector<Gpu::HostVector<int> > host_runtime_int(NumRuntimeIntComps());

            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    Real r[3];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();                
                    p.pos(0) = plo[0] + (iv[0] + r[0])*dx[0];
#if AMREX_SPACEDIM > 1
                    p.pos(1) = plo[1] + (iv[1] + r[1])*dx[1];
#endif
#if AMREX_SPACEDIM > 2
                    p.pos(2) = plo[2] + (iv[2] + r[2])*dx[2];
#endif
                    
                    for (int i = 0; i < NSR; ++i) p.rdata(i) = p.id();
                    for (int i = 0; i < NSI; ++i) p.idata(i) = p.id();
                    
                    host_particles.push_back(p);
                    for (int i = 0; i < NAR; ++i)
                        host_real[i].push_back(p.id());
                    for (int i = 0; i < NAI; ++i)
                        host_int[i].push_back(p.id());
                    for (int i = 0; i < NumRuntimeRealComps(); ++i)
                        host_runtime_real[i].push_back(p.id());
                    for (int i = 0; i < NumRuntimeIntComps(); ++i)
                        host_runtime_int[i].push_back(p.id());
                }
            }
        
            auto& particle_tile = DefineAndReturnParticleTile(lev, mfi.index(), mfi.LocalTileIndex());
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);
            
            Gpu::copy(Gpu::hostToDevice,
                      host_particles.begin(),
                      host_particles.end(),
                      particle_tile.GetArrayOfStructs().begin() + old_size);        

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_real[i].begin(),
                          host_real[i].end(),
                          soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_int[i].begin(),
                          host_int[i].end(),
                          soa.GetIntData(i).begin() + old_size);
            }
            for (int i = 0; i < NumRuntimeRealComps(); ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_runtime_real[i].begin(),
                          host_runtime_real[i].end(),
                          soa.GetRealData(NAR+i).begin() + old_size);
            }

            for (int i = 0; i < NumRuntimeIntComps(); ++i)
            {
                Gpu::copy(Gpu::hostToDevice,
                          host_runtime_int[i].begin(),
                          host_runtime_int[i].end(),
                          soa.GetIntData(NAI+i).begin() + old_size);
            }
        }

        RedistributeLocal();
    }

    void moveParticles (const IntVect& move_dir, int do_random)
    {
        BL_PROFILE("TestParticleContainer::moveParticles");

        for (int lev = 0; lev <= finestLevel(); ++lev)
        {
            const auto dx = Geom(lev).CellSizeArray();
            auto& plev  = GetParticles(lev);
        
            for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                int gid = mfi.index();
                int tid = mfi.LocalTileIndex();            
                auto& ptile = plev[std::make_pair(gid, tid)];
                auto& aos   = ptile.GetArrayOfStructs();
                ParticleType* pstruct = &(aos[0]);            
                const size_t np = aos.numParticles();

                if (do_random == 0)
                {
                    amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
                    {
                        ParticleType& p = pstruct[i];
                        p.pos(0) += move_dir[0]*dx[0];
#if AMREX_SPACEDIM > 1
                        p.pos(1) += move_dir[1]*dx[1];
#endif
#if AMREX_SPACEDIM > 2
                        p.pos(2) += move_dir[2]*dx[2];
#endif
                    });
                }
                else
                {
                    amrex::ParallelFor( np, [=] AMREX_GPU_DEVICE (int i) noexcept
                    {
                        ParticleType& p = pstruct[i];

                        p.pos(0) += (2*amrex::Random()-1)*move_dir[0]*dx[0];
#if AMREX_SPACEDIM > 1
                        p.pos(1) += (2*amrex::Random()-1)*move_dir[1]*dx[1];
#endif
#if AMREX_SPACEDIM > 2
                        p.pos(2) += (2*amrex::Random()-1)*move_dir[2]*dx[2];
#endif
                    });
                }
            }
        }
    }

    void checkAnswer () const
    {
        BL_PROFILE("TestParticleContainer::checkAnswer");
        
        AMREX_ALWAYS_ASSERT(OK());
        
        int num_rr = NumRuntimeRealComps();
        int num_ii = NumRuntimeIntComps();

        for (int lev = 0; lev <= finestLevel(); ++lev)
        {
            auto& plev  = GetParticles(lev);
            for(MFIter mfi = MakeMFIter(lev); mfi.isValid(); ++mfi)
            {
                int gid = mfi.index();
                int tid = mfi.LocalTileIndex();            
                auto& ptile = plev.at(std::make_pair(gid, tid));
                const auto ptd = ptile.getConstParticleTileData();
                const size_t np = ptile.numParticles();
                
                AMREX_FOR_1D ( np, i,
                {
                    for (int j = 0; j < NSR; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_aos[i].rdata(j) == ptd.m_aos[i].id());
                    }
                    for (int j = 0; j < NSI; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_aos[i].idata(j) == ptd.m_aos[i].id());
                    }
                    for (int j = 0; j < NAR; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_rdata[j][i] == ptd.m_aos[i].id());
                    }
                    for (int j = 0; j < NAI; ++j)
                    {
                        AMREX_ALWAYS_ASSERT(ptd.m_idata[j][i] == ptd.m_aos[i].id());
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
    IntVect move_dir;
    int do_random;
    int nsteps;
    int nlevs;
    int do_regrid;
    int sort;
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
    pp.get("move_dir", params.move_dir);
    pp.get("do_random", params.do_random);    
    pp.get("nsteps", params.nsteps);
    pp.get("nlevs", params.nlevs);
    pp.get("do_regrid", params.do_regrid);
    pp.query("num_runtime_real", num_runtime_real);
    pp.query("num_runtime_int", num_runtime_int);

    params.sort = 0;
    pp.query("sort", params.sort);
}

void testRedistribute ()
{
    BL_PROFILE("testRedistribute");
    TestParams params;
    get_test_params(params, "redistribute");

    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = params.is_periodic;

    Vector<IntVect> rr(params.nlevs-1);
    for (int lev = 1; lev < params.nlevs; lev++)
        rr[lev-1] = IntVect(D_DECL(2,2,2));

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
    IntVect lo = IntVect(D_DECL(0, 0, 0));
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
    IntVect nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    amrex::Print() << "About to initialize particles \n";

    pc.InitParticles(nppc);

    pc.checkAnswer();

    auto np_old = pc.TotalNumberOfParticles();

    for (int i = 0; i < params.nsteps; ++i)
    {
        pc.moveParticles(params.move_dir, params.do_random);
        pc.RedistributeLocal();
        if (params.sort) pc.SortParticlesByCell();
        pc.checkAnswer();
    }

    if (params.do_regrid)
    {
        const int NProcs = ParallelDescriptor::NProcs();
        {
            for (int lev = 0; lev < params.nlevs; ++lev)
            {
                DistributionMapping new_dm;
                Vector<int> pmap;
                for (int i = 0; i < ba[lev].size(); ++i) pmap.push_back(i % NProcs);
                new_dm.define(pmap);
                pc.SetParticleDistributionMap(lev, new_dm);
            }
            pc.RedistributeGlobal();
            pc.checkAnswer();
        }

        {
            for (int lev = 0; lev < params.nlevs; ++lev)
            {
                DistributionMapping new_dm;
                Vector<int> pmap;
                for (int i = 0; i < ba[lev].size(); ++i) pmap.push_back((i+1) % NProcs);
                new_dm.define(pmap);
                pc.SetParticleDistributionMap(lev, new_dm);
            }
            pc.RedistributeGlobal();
            pc.checkAnswer();            
        }
    }

    if (geom[0].isAllPeriodic()) AMREX_ALWAYS_ASSERT(np_old == pc.TotalNumberOfParticles());

    // the way this test is set up, if we make it here we pass
    amrex::Print() << "pass \n";
}
