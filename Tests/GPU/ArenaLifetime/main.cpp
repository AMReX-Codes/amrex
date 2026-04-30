#include <AMReX.H>
#include <AMReX_Arena.H>
#include <AMReX_GpuAllocators.H>
#include <AMReX_PODVector.H>

using namespace amrex;

/**
 * Test that PolymorphicArenaAllocator-backed containers can safely outlive
 * amrex::Finalize().
 *
 * This exercises the Arena life-token mechanism: after Finalize destroys the
 * arenas, PolymorphicArenaWrapper falls back to The_Null_Arena (a no-op free)
 * instead of calling free() through a dangling pointer.
 */
int main (int argc, char* argv[])
{
    // A PODVector backed by a polymorphic arena, declared BEFORE Initialize
    // so it outlives Finalize.
    // In more complex situations, e.g., Python workflows, this could equally
    // happen if users create new AMReX containers at runtime and hold on to
    // them, while the original simulation already finalizes earlier than the
    // destruction of this new container occurs.
    PODVector<double, PolymorphicArenaAllocator<double>> vec;

    amrex::Initialize(argc, argv);

    // Point the vector at The_Arena (a CArena on GPU builds, deleted at Finalize)
    vec.setArena(The_Arena());

    // Allocate some data through the arena
    constexpr int N = 1024;
    vec.resize(N, 0.0);

    amrex::Print() << "Allocated " << N << " elements through The_Arena.\n";

    // Finalize destroys The_Arena.  vec still holds pointers into it.
    amrex::Finalize();

    amrex::Print() << "Post-Finalize: vec destructor about to run... ";

    // vec's destructor runs here.  The PolymorphicArenaWrapper should now
    // fall back to The_Null_Arena, whose free() is a no-op.
}
