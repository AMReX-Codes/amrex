#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_Particles.H>

using namespace amrex;

static constexpr int NSR = 0;
static constexpr int NSI = 0;
static constexpr int NAR = 0;
static constexpr int NAI = 2;

void get_position_unit_cell(Real* r, const IntVect& nppc, int i_part)
{
    int nx = nppc[0];
    int ny = nppc[1];
    int nz = nppc[2];
    
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

    using ParticleTileType = ParticleTile<NSR, NSI, NAR, NAI>;

    
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
            std::array<Gpu::HostVector<Real>, NAR> host_real;
            std::array<Gpu::HostVector<int>, NAI> host_int;
            for (IntVect iv = tile_box.smallEnd(); iv <= tile_box.bigEnd(); tile_box.next(iv))
            {
                for (int i_part=0; i_part<num_ppc;i_part++) {
                    Real r[3];
                    get_position_unit_cell(r, a_num_particles_per_cell, i_part);
                
                    Real x = plo[0] + (iv[0] + r[0])*dx[0];
                    Real y = plo[1] + (iv[1] + r[1])*dx[1];
                    Real z = plo[2] + (iv[2] + r[2])*dx[2];
                
                    ParticleType p;
                    p.id()  = ParticleType::NextID();
                    p.cpu() = ParallelDescriptor::MyProc();                
                    p.pos(0) = x;
                    p.pos(1) = y;
                    p.pos(2) = z;
                    
                    for (int i = 0; i < NSR; ++i) p.rdata(i) = i;
                    for (int i = 0; i < NSI; ++i) p.idata(i) = i;
                    
                    host_particles.push_back(p);
                    for (int i = 0; i < NAR; ++i)
                        host_real[i].push_back(i);
                    for (int i = 0; i < NAI; ++i)
                        host_int[i].push_back(i);
                }
            }
        
            auto& particles = GetParticles(lev);
            auto& particle_tile = particles[std::make_pair(mfi.index(), mfi.LocalTileIndex())];
            auto old_size = particle_tile.GetArrayOfStructs().size();
            auto new_size = old_size + host_particles.size();
            particle_tile.resize(new_size);
            
            Gpu::copy(Gpu::hostToDevice, host_particles.begin(), host_particles.end(),
                      particle_tile.GetArrayOfStructs().begin() + old_size);        
            
            auto& soa = particle_tile.GetStructOfArrays();
            for (int i = 0; i < NAR; ++i)
            {
                Gpu::copy(Gpu::hostToDevice, host_real[i].begin(), host_real[i].end(),
                          soa.GetRealData(i).begin() + old_size);
            }
            
            for (int i = 0; i < NAI; ++i)
            {
                Gpu::copy(Gpu::hostToDevice, host_int[i].begin(), host_int[i].end(),
                          soa.GetIntData(i).begin() + old_size);
            }
        }
    }
};

struct Transformer
{    
    int m_factor;

    /**
       \brief This copies the particle but multiplies all the idata by m_factor
       *
     **/
    Transformer (int a_factor)
        : m_factor(a_factor)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst, const SrcData& src,
                     int src_i, int dst_i) const noexcept
    {
        dst.m_aos[dst_i] = src.m_aos[src_i];
        for (int j = 0; j < DstData::NAR; ++j)
            dst.m_rdata[j][dst_i] = src.m_rdata[j][src_i];
        for (int j = 0; j < dst.m_num_runtime_real; ++j)
            dst.m_runtime_rdata[j][dst_i] = src.m_runtime_rdata[j][src_i];
        for (int j = 0; j < DstData::NAI; ++j)
            dst.m_idata[j][dst_i] = m_factor*src.m_idata[j][src_i];
        for (int j = 0; j < dst.m_num_runtime_int; ++j)
            dst.m_runtime_idata[j][dst_i] = src.m_runtime_idata[j][src_i];
    }
};

struct TwoWayTransformer
{    
    int m_factor1;
    int m_factor2;

    /**
       \brief This copies the particle but multiplies all the idata by a factor,
       placing the result in both dst1 and dst2;
       *
     **/
    TwoWayTransformer (int a_factor1, int a_factor2)
        : m_factor1(a_factor1), m_factor2(a_factor2)
    {}

    template <typename DstData, typename SrcData>
    AMREX_GPU_HOST_DEVICE
    void operator() (const DstData& dst1, const DstData& dst2, const SrcData& src,
                     int src_i, int dst1_i, int dst2_i) const noexcept
    {
        dst1.m_aos[dst1_i] = src.m_aos[src_i];
        for (int j = 0; j < DstData::NAR; ++j)
            dst1.m_rdata[j][dst1_i] = src.m_rdata[j][src_i];
        for (int j = 0; j < dst1.m_num_runtime_real; ++j)
            dst1.m_runtime_rdata[j][dst1_i] = src.m_runtime_rdata[j][src_i];
        for (int j = 0; j < DstData::NAI; ++j)
            dst1.m_idata[j][dst1_i] = m_factor1*src.m_idata[j][src_i];
        for (int j = 0; j < dst1.m_num_runtime_int; ++j)
            dst1.m_runtime_idata[j][dst1_i] = src.m_runtime_idata[j][src_i];

        dst2.m_aos[dst2_i] = src.m_aos[src_i];
        for (int j = 0; j < DstData::NAR; ++j)
            dst2.m_rdata[j][dst2_i] = src.m_rdata[j][src_i];
        for (int j = 0; j < dst2.m_num_runtime_real; ++j)
            dst2.m_runtime_rdata[j][dst2_i] = src.m_runtime_rdata[j][src_i];
        for (int j = 0; j < DstData::NAI; ++j)
            dst2.m_idata[j][dst2_i] = m_factor2*src.m_idata[j][src_i];
        for (int j = 0; j < dst2.m_num_runtime_int; ++j)
            dst2.m_runtime_idata[j][dst2_i] = src.m_runtime_idata[j][src_i];
    }
};

struct KeepOddFilter
{    
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        return (src.m_aos[i].id() % 2 == 1);
    }
};

struct KeepEvenFilter
{    
    template <typename SrcData>
    AMREX_GPU_HOST_DEVICE
    int operator() (const SrcData& src, int i) const noexcept
    {
        return (src.m_aos[i].id() % 2 == 0);
    }
};

template <typename PC, typename F>
void transformParticles (PC& pc, F&& f)
{
    BL_PROFILE("transformParticles");

    using ParIter = typename PC::ParConstIterType;
    using ParticleTileType = typename PC::ParticleTileType;
    
    for (int lev = 0; lev <= pc.finestLevel(); ++lev)
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& ptile = pc.ParticlesAt(lev, pti);
            
            ParticleTileType ptile_tmp;
            ptile_tmp.resize(ptile.size());
            
            amrex::transformParticles(ptile_tmp, ptile, std::forward<F>(f));
            ptile.swap(ptile_tmp);
        }
    }
}

template <typename PC, typename F>
void twoWayTransformParticles (PC& dst1, PC& dst2, const PC& src, F&& f)
{
    BL_PROFILE("twoWayTransformParticles");

    using ParIter = typename PC::ParConstIterType;

    for (int lev = 0; lev <= src.finestLevel(); ++lev)
    {
        for(ParIter pti(src, lev); pti.isValid(); ++pti)
        {
            auto& ptile_src = src.ParticlesAt(lev, pti);
            auto& ptile_dst1 = dst1.DefineAndReturnParticleTile(lev, pti);
            auto& ptile_dst2 = dst2.DefineAndReturnParticleTile(lev, pti);
            
            ptile_dst1.resize(ptile_src.size());
            ptile_dst2.resize(ptile_src.size());
            
            amrex::transformParticles(ptile_dst1, ptile_dst2, ptile_src, std::forward<F>(f));
        }
    }
}

template <typename PC>
void testTwoWayTransform (const PC& pc)
{
    using PType = typename PC::SuperParticleType;

    PC pc1(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    pc1.copyParticles(pc);

    auto np_old = pc1.TotalNumberOfParticles();
    
    PC pc2(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    PC pc3(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    
    twoWayTransformParticles(pc2, pc3, pc1, TwoWayTransformer(2, 3));

    auto mx1 = amrex::ReduceMax(pc1, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    auto mx2 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    AMREX_ALWAYS_ASSERT(2*mx1 == mx2);

    auto np_new = pc2.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(np_old == np_new);
    
    auto mx3 = amrex::ReduceMax(pc3, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    AMREX_ALWAYS_ASSERT(mx3 == 3*mx1);

    np_new = pc3.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(np_old == np_new);    
}

template <typename PC, typename F>
void filterParticles (PC& pc, F&& f)
{
    BL_PROFILE("filterParticles");

    using ParIter = typename PC::ParIterType;
    using ParticleTileType = typename PC::ParticleTileType;
    
    for (int lev = 0; lev <= pc.finestLevel(); ++lev)
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& ptile = pc.ParticlesAt(lev, pti);
            
            ParticleTileType ptile_tmp;
            ptile_tmp.resize(ptile.size());
            
            auto num_output = amrex::filterParticles(ptile_tmp, ptile, std::forward<F>(f));

            ptile.swap(ptile_tmp);
            ptile.resize(num_output);
        }
    }
}

template <typename PC, typename Pred, typename F>
void filterAndTransformParticles (PC& pc, Pred&& p, F&& f)
{
    BL_PROFILE("filterAndTransformParticles");

    using ParIter = typename PC::ParIterType;
    using ParticleTileType = typename PC::ParticleTileType;
    
    for (int lev = 0; lev <= pc.finestLevel(); ++lev)
    {
        for(ParIter pti(pc, lev); pti.isValid(); ++pti)
        {
            auto& ptile = pc.ParticlesAt(lev, pti);
            
            ParticleTileType ptile_tmp;
            ptile_tmp.resize(ptile.size());
            
            auto num_output = amrex::filterAndTransformParticles(ptile_tmp, ptile, std::forward<Pred>(p), std::forward<F>(f));
            ptile.swap(ptile_tmp);
            ptile.resize(num_output);
        }
    }
}

template <typename PC, typename Pred, typename F>
void twoWayFilterAndTransformParticles (PC& dst1, PC& dst2, const PC& src, Pred&& p, F&& f)
{
    BL_PROFILE("twoWayFilterAndTransformParticles");

    using ParIter = typename PC::ParConstIterType;    
    for (int lev = 0; lev <= src.finestLevel(); ++lev)
    {
        for(ParIter pti(src, lev); pti.isValid(); ++pti)
        {
            auto& ptile_src = src.ParticlesAt(lev, pti);
            auto& ptile_dst1 = dst1.DefineAndReturnParticleTile(lev, pti);
            auto& ptile_dst2 = dst2.DefineAndReturnParticleTile(lev, pti);
            
            ptile_dst1.resize(ptile_src.size());
            ptile_dst2.resize(ptile_src.size());
            
            auto num_output = amrex::filterAndTransformParticles(ptile_dst1,
                                                                 ptile_dst2,
                                                                 ptile_src,
                                                                 std::forward<Pred>(p),
                                                                 std::forward<F>(f));
            ptile_dst1.resize(num_output);
            ptile_dst2.resize(num_output);
        }
    }
}

template <typename PC>
void testTransform (const PC& pc)
{
    using PType = typename PC::SuperParticleType;

    PC pc2(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    pc2.copyParticles(pc);

    auto np_old = pc2.TotalNumberOfParticles();
    
    auto mx1 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    transformParticles(pc2, Transformer(2));
    
    auto mx2 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    AMREX_ALWAYS_ASSERT(2*mx1 == mx2);

    auto np_new = pc2.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(np_old == np_new);
    
    pc2.clearParticles();
    pc2.copyParticles(pc);
    transformParticles(pc2, Transformer(3));
    
    auto mx3 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    AMREX_ALWAYS_ASSERT(mx3 == 3*mx1);    
}

template <typename PC>
void testFilter (const PC& pc)
{
    using PType = typename PC::SuperParticleType;

    auto np_old = pc.TotalNumberOfParticles();

    PC pc2(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    pc2.copyParticles(pc, KeepOddFilter());

    auto np_new = pc2.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(2*np_new == np_old);

    auto all_odd = amrex::ReduceLogicalAnd(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.id() % 2 == 1; });

    AMREX_ALWAYS_ASSERT(all_odd);

    pc2.clearParticles();
    pc2.copyParticles(pc, KeepEvenFilter());

    np_new = pc2.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(2*np_new == np_old);

    auto all_even = amrex::ReduceLogicalAnd(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.id() % 2 == 0; });

    AMREX_ALWAYS_ASSERT(all_even);

    filterParticles(pc2, KeepOddFilter());

    np_new = pc2.TotalNumberOfParticles();

    AMREX_ALWAYS_ASSERT(np_new == 0);
}

template <typename PC>
void testFilterAndTransform (const PC& pc)
{
    using PType = typename PC::SuperParticleType;

    PC pc2(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    pc2.copyParticles(pc);

    auto mx1 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    auto np_old = pc2.TotalNumberOfParticles();
    
    filterAndTransformParticles(pc2, KeepOddFilter(), Transformer(3));

    auto mx2 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });
    
    auto np_new = pc2.TotalNumberOfParticles();
    
    AMREX_ALWAYS_ASSERT(2*np_new == np_old);
    AMREX_ALWAYS_ASSERT(mx2 == 3*mx1);
}

template <typename PC>
void testTwoWayFilterAndTransform (const PC& pc)
{
    using PType = typename PC::SuperParticleType;

    PC pc1(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    pc1.copyParticles(pc);

    auto np_old = pc1.TotalNumberOfParticles();
    
    PC pc2(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    PC pc3(pc.Geom(0), pc.ParticleDistributionMap(0), pc.ParticleBoxArray(0));
    
    twoWayFilterAndTransformParticles(pc2, pc3, pc1, KeepEvenFilter(), TwoWayTransformer(2, 3));
    
    auto mx1 = amrex::ReduceMax(pc1, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    auto mx2 = amrex::ReduceMax(pc2, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });

    auto mx3 = amrex::ReduceMax(pc3, [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> int { return p.idata(NSI+1); });
    
    auto np_new2 = pc2.TotalNumberOfParticles();
    auto np_new3 = pc2.TotalNumberOfParticles();
    
    AMREX_ALWAYS_ASSERT(2*np_new2 == np_old);
    AMREX_ALWAYS_ASSERT(2*np_new3 == np_old);
    AMREX_ALWAYS_ASSERT(mx2 == 2*mx1);
    AMREX_ALWAYS_ASSERT(mx3 == 3*mx1);
}

struct TestParams
{
    IntVect size;
    int max_grid_size;
    int num_ppc;
};

void testTransformations ();

int main (int argc, char* argv[])
{
    amrex::Initialize(argc,argv);

    amrex::Print() << "Running particle reduction test \n";
    testTransformations();

    amrex::Finalize();
}

void get_test_params (TestParams& params, const std::string& prefix)
{
    ParmParse pp(prefix);
    pp.get("size", params.size);
    pp.get("max_grid_size", params.max_grid_size);
    pp.get("num_ppc", params.num_ppc);
}

void testTransformations ()
{
    BL_PROFILE("testTransformations");
    TestParams params;
    get_test_params(params, "reduce");

    RealBox real_box;
    for (int n = 0; n < BL_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, params.size[n]);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(params.size[0]-1,params.size[1]-1,params.size[2]-1));
    const Box domain(domain_lo, domain_hi);

    int coord = 0;
    int is_per[BL_SPACEDIM];
    for (int i = 0; i < BL_SPACEDIM; i++)
        is_per[i] = 1;
    Geometry geom(domain, &real_box, coord, is_per);
    
    BoxArray ba(domain);
    ba.maxSize(params.max_grid_size);
    DistributionMapping dm(ba);

    TestParticleContainer pc(geom, dm, ba);

    int npc = params.num_ppc;
    IntVect nppc = IntVect(AMREX_D_DECL(npc, npc, npc));

    if (ParallelDescriptor::MyProc() == dm[0])
        amrex::Print() << "About to initialize particles \n";

    pc.InitParticles(nppc);

    testTransform(pc);

    testFilter(pc);

    testFilterAndTransform(pc);

    testTwoWayTransform(pc);

    testTwoWayFilterAndTransform(pc);
    
    amrex::Print() << "pass \n";
}
