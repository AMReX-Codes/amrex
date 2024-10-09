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

void addParticles ()
{
    using PC = ParticleContainerPureSoA<AMREX_SPACEDIM+1, 2>;
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

    PC pc(geom, dm, ba);

    amrex::Print() << "Original Real SoA component names are: ";
    for (auto& n : pc.GetRealSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    amrex::Print() << "Original Int SoA component names are: ";
    for (auto& n : pc.GetIntSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    amrex::Print() << "Adding runtime comps. \n";
    pc.AddRealComp("real_comp1");
    pc.AddRealComp(); // without name - should be real_comp2
    pc.AddIntComp(); // without name - should be int_comp0

    amrex::Print() << "New Real SoA component names are: ";
    for (auto& n : pc.GetRealSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    amrex::Print() << "New Int SoA component names are: ";
    for (auto& n : pc.GetIntSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    amrex::Print() << "Reset compile-time SoA names \n";
    pc.SetSoACompileTimeNames({AMREX_D_DECL("x", "y", "z"), "w"}, {"i1", "i2"});

    amrex::Print() << "New Real SoA component names are: ";
    for (auto& n : pc.GetRealSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    amrex::Print() << "New Int SoA component names are: ";
    for (auto& n : pc.GetIntSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";

    int const NArrayReal = PC::NArrayReal;
    int const NArrayInt = PC::NArrayInt;
    using ParticleType = typename PC::ParticleType;

    const int add_num_particles = 5;
    auto& ptile1 = pc.DefineAndReturnParticleTile(0, 0, 0);
    ptile1.resize(add_num_particles);

    for (int i = 0; i < add_num_particles; ++i)
    {
        for (int d = 0; d < AMREX_SPACEDIM; d++) {
            ptile1.pos(i, d) = 12.0;
        }
        ptile1.getParticleTileData().rdata(AMREX_SPACEDIM)[i] = 1.2;  // w

        ptile1.push_back_int(0, int(ParticleType::NextID()));
        ptile1.push_back_int(1, amrex::ParallelDescriptor::MyProc());
    }

    int lev=0;
    using MyParIter = ParIter_impl<ParticleType, NArrayReal, NArrayInt>;
    for (MyParIter pti(pc, lev); pti.isValid(); ++pti) {
        auto& soa = pti.GetStructOfArrays();
        AMREX_D_TERM(
                     auto *xp = soa.GetRealData("x").data();,
                     auto *yp = soa.GetRealData("y").data();,
                     auto *zp = soa.GetRealData("z").data();
                     );
        auto *wp = soa.GetRealData("w").data();

        const int np = pti.numParticles();
        ParallelFor( np, [=] AMREX_GPU_DEVICE (long ip)
        {
            AMREX_D_TERM(
                         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(xp[ip] == 12_prt,
                                                          "pos attribute expected to be 12");,
                         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(yp[ip] == 12_prt,
                                                          "pos attribute expected to be 12");,
                         AMREX_ALWAYS_ASSERT_WITH_MESSAGE(zp[ip] == 12_prt,
                                                          "pos attribute expected to be 12");
                         );
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(wp[ip] == 1.2_prt,
                                             "pos attribute expected to be 1.2");
        });
    }
}

int main (int argc, char* argv[])
 {
    amrex::Initialize(argc,argv);
    {
        addParticles();
    }
    amrex::Finalize();
 }
