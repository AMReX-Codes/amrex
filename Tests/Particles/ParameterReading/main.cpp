#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>

using namespace amrex;

// Two particle container types with different template arguments.
using ContainerA = ParticleContainerPureSoA<AMREX_SPACEDIM, 0>;
using ContainerB = ParticleContainer<1, 0>;

namespace {

Geometry make_geometry ()
{
    RealBox rb(AMREX_D_DECL(0.0, 0.0, 0.0), AMREX_D_DECL(1.0, 1.0, 1.0));
    Box domain(IntVect(0), IntVect(15));
    return Geometry(domain, rb, CoordSys::cartesian, {AMREX_D_DECL(0, 0, 0)});
}

void check_defaults ()
{
    AMREX_ALWAYS_ASSERT(!ParticleContainerBase::do_tiling);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::memEfficientSort);
    AMREX_ALWAYS_ASSERT(!ParticleContainerBase::use_comms_arena);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::redistribute_early_out);
}

// First cycle: no particles.* inputs, so all parameters start at their defaults.
void values_set_in_code_are_kept ()
{
    Geometry geom = make_geometry();
    BoxArray ba(geom.Domain());
    DistributionMapping dm(ba);

    ContainerA pa(geom, dm, ba);
    check_defaults();

    // Change the process-wide parameters in code after the first container is built ...
    ParticleContainerBase::do_tiling = true;
    ParticleContainerBase::tile_size = IntVect(AMREX_D_DECL(4, 4, 4));
    ParticleContainerBase::memEfficientSort = false;
    ParticleContainerBase::use_comms_arena = true;
    ParticleContainerBase::redistribute_early_out = false;

    // ... and building a container of another type must not undo that.
    ContainerB pb(geom, dm, ba);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::do_tiling);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::tile_size == IntVect(AMREX_D_DECL(4, 4, 4)));
    AMREX_ALWAYS_ASSERT(!ParticleContainerBase::memEfficientSort);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::use_comms_arena);
    AMREX_ALWAYS_ASSERT(!ParticleContainerBase::redistribute_early_out);
}

// Second cycle: particles.* inputs are set, and Finalize has restored the defaults.
void inputs_are_read_after_reinitialize ()
{
    // Finalize restored the defaults; nothing has been read yet in this cycle.
    check_defaults();

    Geometry geom = make_geometry();
    BoxArray ba(geom.Domain());
    DistributionMapping dm(ba);

    ContainerA pa(geom, dm, ba);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::do_tiling);
    AMREX_ALWAYS_ASSERT(!ParticleContainerBase::memEfficientSort);
    AMREX_ALWAYS_ASSERT(ParticleContainerBase::redistribute_early_out);
    AMREX_ALWAYS_ASSERT(pa.GetUsePrePost());
    AMREX_ALWAYS_ASSERT(!pa.GetUseUnlink());

    // Per-container parameters are read for every container, not just the first of a type.
    ContainerA pa2(geom, dm, ba);
    AMREX_ALWAYS_ASSERT(pa2.GetUsePrePost());
    AMREX_ALWAYS_ASSERT(!pa2.GetUseUnlink());
}

}

int main (int argc, char* argv[])
{
    // MPI outlives both AMReX cycles, so AMReX must not finalize it.
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    amrex::Initialize(argc, argv);
    values_set_in_code_are_kept();
    amrex::Finalize();

    amrex::Initialize(argc, argv, true, MPI_COMM_WORLD, [] () {
        ParmParse pp("particles");
        pp.add("do_tiling", true);
        pp.add("do_mem_efficient_sort", false);
        pp.add("use_prepost", true);
        pp.add("do_unlink", false);
    });
    inputs_are_read_after_reinitialize();
    amrex::Print() << "ParameterReading test passed\n";
    amrex::Finalize();

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
