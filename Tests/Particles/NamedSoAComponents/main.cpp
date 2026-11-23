#include <AMReX.H>
#include <AMReX_Particle.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_ParticleTile.H>
#include <AMReX_ParIter.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_GpuContainers.H>

#include <array>
#include <stdexcept>

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
        ptile1.getParticleTileData().rdata(AMREX_SPACEDIM)[i] = ParticleReal(1.2);  // w

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

    // The runtime component names must stay in sync with the number of
    // components when the runtime components are resized.
    amrex::Print() << "Resize runtime comps. \n";
    auto check_names = [&] () {
        AMREX_ALWAYS_ASSERT(int(pc.GetRealSoANames().size()) == pc.NumRealComps());
        AMREX_ALWAYS_ASSERT(int(pc.GetIntSoANames().size()) == pc.NumIntComps());
    };

    // an empty tile, which ParIter does not visit, must be resized too
    auto& empty_tile = pc.DefineAndReturnParticleTile(0, 1, 0);
    AMREX_ALWAYS_ASSERT(empty_tile.numParticles() == 0);
    auto check_empty_tile = [&] () {
        AMREX_ALWAYS_ASSERT(empty_tile.NumRealComps() == pc.NumRealComps());
        AMREX_ALWAYS_ASSERT(empty_tile.NumIntComps() == pc.NumIntComps());
    };

    // grow: new components get default names that can be looked up
    pc.ResizeRuntimeRealComp(4, true);
    pc.ResizeRuntimeIntComp(3, true);
    check_names();
    check_empty_tile();
    AMREX_ALWAYS_ASSERT(pc.GetRealCompIndex("real_comp4") == NArrayReal + 3);
    AMREX_ALWAYS_ASSERT(pc.GetIntCompIndex("int_comp4") == NArrayInt + 2);
    for (MyParIter pti(pc, lev); pti.isValid(); ++pti) {
        auto& soa = pti.GetStructOfArrays();
        const auto np = static_cast<std::size_t>(pti.numParticles());
        AMREX_ALWAYS_ASSERT(soa.GetRealData("real_comp4").size() == np);
        AMREX_ALWAYS_ASSERT(soa.GetIntData("int_comp4").size() == np);
    }

    // shrink: names of removed components are dropped ...
    pc.ResizeRuntimeRealComp(1, true);
    pc.ResizeRuntimeIntComp(0, true);
    check_names();
    check_empty_tile();
    AMREX_ALWAYS_ASSERT(pc.HasRealComp("real_comp1"));
    AMREX_ALWAYS_ASSERT(!pc.HasRealComp("real_comp2"));
    AMREX_ALWAYS_ASSERT(!pc.HasIntComp("int_comp2"));

    // ... so adding a component with a default name works again
    pc.AddRealComp(); // real_comp2
    pc.AddIntComp();  // int_comp2
    check_names();
    AMREX_ALWAYS_ASSERT(pc.GetRealCompIndex("real_comp2") == NArrayReal + 1);
    AMREX_ALWAYS_ASSERT(pc.GetIntCompIndex("int_comp2") == NArrayInt);

    // grow: a default name that is already in use is an error, and the
    // container is left unchanged
    pc.AddRealComp("real_comp4"); // takes the slot whose default name is real_comp3
    auto const names_before = pc.GetRealSoANames();
    bool threw = false;
    try {
        pc.ResizeRuntimeRealComp(4, true); // would add a second real_comp4
    } catch (std::runtime_error const&) {
        threw = true;
    }
    AMREX_ALWAYS_ASSERT(threw);
    AMREX_ALWAYS_ASSERT(pc.NumRuntimeRealComps() == 3);
    AMREX_ALWAYS_ASSERT(pc.GetRealSoANames() == names_before);
    check_names();
    check_empty_tile();

    amrex::Print() << "Final Real SoA component names are: ";
    for (auto& n : pc.GetRealSoANames()) {
        amrex::Print() << n << ", ";
    }
    amrex::Print() << "\n";
}

int main (int argc, char* argv[])
 {
    {
        amrex::ParmParse pp("amrex");
        pp.add("the_arena_is_managed", 1);
    }
    amrex::Initialize(argc,argv);
    {
        addParticles();
    }
    amrex::Finalize();
 }
