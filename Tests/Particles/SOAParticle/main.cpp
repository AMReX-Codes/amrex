#include <AMReX.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParIter.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <array>

using namespace amrex;

template <typename T_PC,template<class> class Allocator=DefaultAllocator>
void addParticles ()
{
    int is_per[AMREX_SPACEDIM];
    for (int & d : is_per) {
        d = 1;
    }

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++)
    {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 100.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(127, 127, 127));
    const Box base_domain(domain_lo, domain_hi);

    Geometry geom(base_domain, &real_box, CoordSys::cartesian, is_per);
    BoxArray ba(base_domain);
    ba.maxSize(64);

    DistributionMapping dm(ba);

    T_PC pc(geom, dm, ba);
    int const NArrayReal = pc.NArrayReal;
    int const NArrayInt = pc.NArrayInt;

    using ParticleType = typename T_PC::ParticleType;
    using ParticleTileDataType = typename T_PC::ParticleTileType::ParticleTileDataType;

    const int add_num_particles = 5;

    auto& ptile1 = pc.DefineAndReturnParticleTile(0, 0, 0);
    ptile1.resize(add_num_particles);

    for (int i = 0; i < add_num_particles; ++i)
    {
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            ptile1.pos(i, d) = 12.0;
        }
        ptile1.getParticleTileData().rdata(AMREX_SPACEDIM)[i] = 1.2;  // w

        ptile1.push_back_int(0, ParticleType::NextID());
        ptile1.push_back_int(1, amrex::ParallelDescriptor::MyProc());
    }

    int lev=0;
    // int numparticles=0;
    using MyParIter = ParIter_impl<ParticleType, NArrayReal, NArrayInt>;
    for (MyParIter pti(pc, lev); pti.isValid(); ++pti) {
        const int np = pti.numParticles();
        // preparing access to particle data: SoA of Reals
        auto& soa = pti.GetStructOfArrays();
        auto soa_real = soa.GetRealData();
        auto size = soa.size();
        amrex::ParticleReal* const AMREX_RESTRICT part_x = soa_real[0].dataPtr();
        amrex::ParticleReal* const AMREX_RESTRICT part_y = AMREX_SPACEDIM >= 2 ? soa_real[1].dataPtr() : nullptr;
        amrex::ParticleReal* const AMREX_RESTRICT part_z = AMREX_SPACEDIM >= 3 ? soa_real[2].dataPtr() : nullptr;
        amrex::ParticleReal* const AMREX_RESTRICT part_w = soa_real[AMREX_SPACEDIM].dataPtr();
        auto& soa_int = pti.GetStructOfArrays().GetIntData();
        amrex::ignore_unused(size, part_x, part_y, part_z, part_w, soa_int);

        // Iterating over old Particles
        // ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        // {
        //     ParticleType& AMREX_RESTRICT p = aos_ptr[ip];
        //     p.pos(0) += 1;
        //     p.pos(1) += 1;
        //     p.pos(2) += 1;

        //     amrex::ParticleReal & AMREX_RESTRICT x = part_x[ip];
        //     amrex::ParticleReal & AMREX_RESTRICT y = part_y[ip];
        //     amrex::ParticleReal & AMREX_RESTRICT z = part_z[ip];
        //     amrex::ParticleReal & AMREX_RESTRICT a = part_aaa[ip];

        //     x += 1.0;
        //     y += 1.0;
        //     z += 1.0;
        //     a += 1.0;
        // });

        // Iterating over SoA Particles
        ParticleTileDataType ptd = pti.GetParticleTile().getParticleTileData();

        ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            ParticleType p(ptd, ip);
            for (int d = 0; d < AMREX_SPACEDIM; d++) {
                p.pos(d) += 1_prt;
                AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ptd.rdata(d)[ip] == 13_prt,
                                                 "pos attribute expected to be 13");
            }

            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ptd.rdata(AMREX_SPACEDIM)[ip] == 1.2_prt,
                                             "w attribute expected to be 1.2");
        });


    }

    // create a host-side particle buffer
    auto tmp = pc.template make_alike<amrex::PinnedArenaAllocator>();
    tmp.copyParticles(pc, true);

    using MyPinnedParIter = ParIter_impl<ParticleType, NArrayReal, NArrayInt, amrex::PinnedArenaAllocator>;

    for (MyPinnedParIter pti(tmp, lev); pti.isValid(); ++pti) {
        auto& particle_attributes = pti.GetStructOfArrays();
        auto& real_comp0 = particle_attributes.GetRealData(0);
        auto&  int_comp1  = particle_attributes.GetIntData(1);
        for (int i = 0; i < pti.numParticles(); ++i) {
            real_comp0[i] += 1;
            int_comp1[i] += 1;
        }
    }

    tmp.Redistribute();

    using ConstPTDType = typename T_PC::ParticleTileType::ConstParticleTileDataType;
    amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
    auto r = amrex::ParticleReduce<
        amrex::ReduceData<
            amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal, amrex::ParticleReal, amrex::ParticleReal,
            amrex::ParticleReal>
    >(
        pc,
        [=] AMREX_GPU_DEVICE(const ConstPTDType& ptd, const int i) noexcept
        {
            amrex::ParticleReal const x = ptd.rdata(0)[i];
            amrex::ParticleReal const y = AMREX_SPACEDIM >= 2 ? ptd.rdata(1)[i] : 0.0;
            amrex::ParticleReal const z = AMREX_SPACEDIM >= 3 ? ptd.rdata(2)[i] : 0.0;

            amrex::ParticleReal const w = ptd.rdata(AMREX_SPACEDIM)[i];

            return amrex::makeTuple(x, x*x, y, y*y, z, z*z, w);
        },
        reduce_ops
    );
    amrex::ignore_unused(r);

    // Reduce for SoA Particle Struct
    /*
    using PTDType = typename T_PC::ParticleTileType::ConstParticleTileDataType;
    amrex::ReduceOps<ReduceOpSum, ReduceOpMin, ReduceOpMax> reduce_ops;
    auto r = amrex::ParticleReduce<ReduceData<amrex::Real, amrex::Real,int>> (
                 pc, [=] AMREX_GPU_DEVICE (const PTDType& ptd, const int i) noexcept
                               -> amrex::GpuTuple<amrex::Real,amrex::Real,int>
             {

                const amrex::Real a = ptd.rdata(1)[i];
                const amrex::Real b = ptd.rdata(2)[i];
                const int c = ptd.idata(1)[i];
                return {a, b, c};
             }, reduce_ops);

    AMREX_ALWAYS_ASSERT(amrex::get<0>(r) == amrex::Real(std::pow(256, AMREX_SPACEDIM)));
    AMREX_ALWAYS_ASSERT(amrex::get<1>(r) == 2.0);
    AMREX_ALWAYS_ASSERT(amrex::get<2>(r) == 1);
    */
}




int main(int argc, char* argv[])
 {
    amrex::Initialize(argc,argv);
    {
        addParticles< ParticleContainerPureSoA<4, 2> > ();
    }
    amrex::Finalize();
 }
