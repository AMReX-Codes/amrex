#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

#include <cmath>

using namespace amrex;

static constexpr int NSR = 4;
static constexpr int NSI = 3;
static constexpr int NAR = 2;
static constexpr int NAI = 1;

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
{
        int nx = nppc[0];
#if AMREX_SPACEDIM >= 2
        int ny = nppc[1];
#else
        int ny = 1;
#endif
#if AMREX_SPACEDIM == 3
        int nz = nppc[2];
#else
        int nz = 1;
#endif

        AMREX_D_TERM(int ix_part = i_part/(ny * nz);,
                     int iy_part = (i_part % (ny * nz)) % ny;,
                     int iz_part = (i_part % (ny * nz)) / ny;)

        AMREX_D_TERM(r[0] = (Real(0.5)+Real(ix_part))/Real(nx);,
                     r[1] = (Real(0.5)+Real(iy_part))/Real(ny);,
                     r[2] = (Real(0.5)+Real(iz_part))/Real(nz);)
}

class TestParticleContainer
    : public amrex::ParticleContainer<NSR, NSI, NAR, NAI>
{

public:

    TestParticleContainer (const amrex::Geometry            & a_geom,
                           const amrex::DistributionMapping & a_dmap,
                           const amrex::BoxArray            & a_ba)
        : amrex::ParticleContainer<NSR, NSI, NAR, NAI>(a_geom, a_dmap, a_ba)
    {}

    void InitParticles (const amrex::IntVect& a_num_particles_per_cell)
    {
        BL_PROFILE("InitParticles");
        const int lev = 0;
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
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    Real r[AMREX_SPACEDIM];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);

                    AMREX_D_TERM(auto x = static_cast<ParticleReal> (plo[0] + (iv[0] + r[0])*dx[0]);,
                                 auto y = static_cast<ParticleReal> (plo[1] + (iv[1] + r[1])*dx[1]);,
                                 auto z = static_cast<ParticleReal> (plo[2] + (iv[2] + r[2])*dx[2]);)

                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();
                    AMREX_D_TERM(p.pos(0) = x;,
                                 p.pos(1) = y;,
                                 p.pos(2) = z;)

                    for (int i = 0; i < NSR; ++i) { p.rdata(i) = ParticleReal(i); }
                    for (int i = 0; i < NSI; ++i) { p.idata(i) = i; }

                    host_particles.push_back(p);
                    for (int i = 0; i < NAR; ++i) {
                        host_real[i].push_back(ParticleReal(i));
                    }
                    for (int i = 0; i < NAI; ++i) {
                        host_int[i].push_back(i);
                    }
                }
            }

            auto& particles = GetParticles(lev);
            auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);

            Gpu::copyAsync(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                           particle_tile.GetArrayOfStructs().begin() + old_size);

            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice, host_real[i].begin(), host_real[i].end(),
                               soa.GetRealData(i).begin() + old_size);
            }

            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copyAsync(Gpu::hostToDevice, host_int[i].begin(), host_int[i].end(),
                               soa.GetIntData(i).begin() + old_size);
            }

            Gpu::streamSynchronize();
        }
    }
};

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int num_ppc;
};

void testReduce();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running particle reduction test \n";
    testReduce();

    amrex::Finalize();
}

void get_test_params(TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
}

void testReduce ()
{
    BL_PROFILE("testReduce");
    TestParams params;
    get_test_params(params, "reduce");

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, Real(0.0));
        real_box.setHi(n, Real(params.size[n]));
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[] = {AMREX_D_DECL(1,1,1)};
    Geometry geom(domain, &real_box, coord, is_per);

    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    TestParticleContainer pc(geom, dm, ba);

    IntVect nppc(params.num_ppc);

    if (ParallelDescriptor::MyProc() == dm[0]) {
        amrex::Print() << "About to initialize particles \n";
    }

    pc.InitParticles(nppc);

    using SPType  = typename TestParticleContainer::SuperParticleType;
    using PType   = typename TestParticleContainer::ParticleType;
    using PTDType = typename TestParticleContainer::ParticleTileType::ConstParticleTileDataType;

    auto sm = amrex::ReduceSum(pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> ParticleReal { return p.rdata(1); });
    AMREX_ALWAYS_ASSERT(sm == pc.TotalNumberOfParticles());

    auto sm2 = amrex::ReduceSum(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> ParticleReal { return -p.rdata(NSR+1); });
    AMREX_ALWAYS_ASSERT(sm2 == -pc.TotalNumberOfParticles());

    auto mn = amrex::ReduceMin(pc, [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> ParticleReal { return ptd.m_aos[i].rdata(1);});
    AMREX_ALWAYS_ASSERT(mn == 1);

    auto mn2 = amrex::ReduceMin(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> ParticleReal { return p.rdata(NSR+1); });
    AMREX_ALWAYS_ASSERT(mn2 == 1);

    auto mx = amrex::ReduceMax(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> ParticleReal { return p.rdata(1); });
    AMREX_ALWAYS_ASSERT(mx == 1);

    auto mx2 = amrex::ReduceMax(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> int { return p.idata(NSI); });
    AMREX_ALWAYS_ASSERT(mx2 == 0);

    {
        auto r = amrex::ReduceLogicalOr(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> int { return p.id() == 1; });
        AMREX_ALWAYS_ASSERT(r == 1);
    }

    {
        auto r = amrex::ReduceLogicalOr(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> int { return p.id() == -1; });
        AMREX_ALWAYS_ASSERT(r == 0);
    }

    {
        auto r = amrex::ReduceLogicalAnd(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> int { return p.id() == p.id(); });
        AMREX_ALWAYS_ASSERT(r == 1);
    }

    {
        auto r = amrex::ReduceLogicalAnd(pc, [=] AMREX_GPU_HOST_DEVICE (const SPType& p) -> int { return p.id() == 1; });
        AMREX_ALWAYS_ASSERT(r == 0);
    }

    {
        amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
        auto r = amrex::ParticleReduce<ReduceData<amrex::ParticleReal, amrex::ParticleReal,int>> (
         pc, [=] AMREX_GPU_DEVICE (const SPType& p) noexcept -> amrex::GpuTuple<amrex::ParticleReal,amrex::ParticleReal,int>
           {
               const amrex::ParticleReal a = p.rdata(1);
               const amrex::ParticleReal b = p.rdata(2);
               const int c = p.idata(1);
               return {a, b, c};
           }, reduce_ops);

        AMREX_ALWAYS_ASSERT(amrex::get<0>(r) == amrex::Real(std::pow(256, AMREX_SPACEDIM)));
        AMREX_ALWAYS_ASSERT(amrex::get<1>(r) == 2.0);
        AMREX_ALWAYS_ASSERT(amrex::get<2>(r) == 1);
    }

    {
        // Tiles whose particle counts are not multiples of 8, so that the CPU
        // reductions also go through their remainder loop.
        const Box small_domain(IntVect(0), IntVect(4));
        RealBox small_real_box({AMREX_D_DECL(0., 0., 0.)}, {AMREX_D_DECL(5., 5., 5.)});
        Geometry small_geom(small_domain, &small_real_box, coord, is_per);
        BoxArray small_ba(small_domain);
        small_ba.maxSize(3);
        DistributionMapping small_dm(small_ba);

        TestParticleContainer small_pc(small_geom, small_dm, small_ba);
        small_pc.InitParticles(IntVect(1));
        const auto np = static_cast<ParticleReal>(small_pc.TotalNumberOfParticles());

        // one particle per cell at x = i + 0.5 for i = 0, ..., 4
        auto xsum = amrex::ReduceSum(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> ParticleReal { return p.pos(0); });
        AMREX_ALWAYS_ASSERT(xsum == ParticleReal(2.5)*np);

        auto xmin = amrex::ReduceMin(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> ParticleReal { return p.pos(0); });
        AMREX_ALWAYS_ASSERT(xmin == ParticleReal(0.5));

        auto xmax = amrex::ReduceMax(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> ParticleReal { return p.pos(0); });
        AMREX_ALWAYS_ASSERT(xmax == ParticleReal(4.5));

        auto any_hi = amrex::ReduceLogicalOr(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.pos(0) > ParticleReal(4.); });
        AMREX_ALWAYS_ASSERT(any_hi == 1);

        auto all_hi = amrex::ReduceLogicalAnd(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.pos(0) > ParticleReal(4.); });
        AMREX_ALWAYS_ASSERT(all_hi == 0);

        auto all_in = amrex::ReduceLogicalAnd(small_pc, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.pos(0) > ParticleReal(0.); });
        AMREX_ALWAYS_ASSERT(all_in == 1);
    }

    amrex::Print() << "pass \n";
}
