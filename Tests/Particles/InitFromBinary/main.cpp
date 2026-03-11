#include <AMReX.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Particles.H>

using namespace amrex;

void test_init_binary ()
{
    // The binary file 32.nyx has NP=32768, DM=3, NX=4
    // with particle positions in [0, 8]^3
    constexpr int NStructReal = 4;
    constexpr int NStructInt  = 0;
    constexpr int NArrayReal  = 0;
    constexpr int NArrayInt   = 0;

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; n++) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 8.0);
    }

    IntVect domain_lo(AMREX_D_DECL(0, 0, 0));
    IntVect domain_hi(AMREX_D_DECL(31, 31, 31));
    const Box domain(domain_lo, domain_hi);

    int is_per[] = {AMREX_D_DECL(1, 1, 1)};
    Geometry geom(domain, &real_box, CoordSys::cartesian, is_per);

    BoxArray ba(domain);
    ba.maxSize(16);

    DistributionMapping dmap(ba);

    using MyPC = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;
    MyPC pc(geom, dmap, ba);

    pc.InitFromBinaryFile("32.nyx", NStructReal);

    Long np = pc.TotalNumberOfParticles();
    amrex::Print() << "Total number of particles: " << np << "\n";
    AMREX_ALWAYS_ASSERT(np == 32768);

    // Particle masses in this file are 1.0, 2.0, ..., 32768.0
    // so the total mass should be 32768 * 32769 / 2 = 536887296
    using PType = MyPC::SuperParticleType;
    Real total_mass = amrex::ReduceSum(pc,
                                       [=] AMREX_GPU_HOST_DEVICE (const PType& p) -> Real
                                       {
                                           return p.rdata(0);
                                       });

    ParallelDescriptor::ReduceRealSum(total_mass);

    amrex::Print() << "Total particle mass: " << total_mass << "\n";
    AMREX_ALWAYS_ASSERT(total_mass == 536887296.0);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    test_init_binary();
    amrex::Finalize();
}
