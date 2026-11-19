/**
 * Test for the HDF5 particle checkpoint "pre/post" path.
 *
 * Writes the same particle data twice, once with the regular
 * CheckpointHDF5() call and once wrapped in CheckpointPreHDF5() /
 * CheckpointPostHDF5() with SetUsePrePost(true). Verifies that
 *   1. both HDF5 files have the same size, i.e. CheckpointPostHDF5()
 *      does not append anything to the HDF5 file, and
 *   2. restarting from the pre/post checkpoint recovers the original
 *      particles.
 */
#include <AMReX.H>
#include <AMReX_Particles.H>

#include <filesystem>

using namespace amrex;

constexpr int NStructReal = 2;
constexpr int NStructInt  = 1;
constexpr int NArrayReal  = 2;
constexpr int NArrayInt   = 1;

using MyPC = ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>;

void verify_same (MyPC& pc_orig, MyPC& pc_new)
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        pc_orig.TotalNumberOfParticles() == pc_new.TotalNumberOfParticles(),
        "Particle count mismatch after restart");

    using PTDType = typename MyPC::ParticleTileType::ConstParticleTileDataType;

    for (int icomp = 0; icomp < NStructReal + NArrayReal; ++icomp) {
        auto f = [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Real {
            return (icomp < NStructReal) ? Real(ptd.m_aos[i].rdata(icomp))
                                         : Real(ptd.m_rdata[icomp-NStructReal][i]);
        };
        auto sm_orig = amrex::ReduceSum(pc_orig, f);
        auto sm_new  = amrex::ReduceSum(pc_new, f);
        ParallelDescriptor::ReduceRealSum(sm_orig);
        ParallelDescriptor::ReduceRealSum(sm_new);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(amrex::almostEqual(sm_orig, sm_new),
            "Real component sum mismatch after restart (comp " + std::to_string(icomp) + ")");
    }

    for (int icomp = 0; icomp < NStructInt + NArrayInt; ++icomp) {
        auto f = [=] AMREX_GPU_HOST_DEVICE (const PTDType& ptd, const int i) -> Long {
            return (icomp < NStructInt) ? Long(ptd.m_aos[i].idata(icomp))
                                        : Long(ptd.m_idata[icomp-NStructInt][i]);
        };
        auto sm_orig = amrex::ReduceSum(pc_orig, f);
        auto sm_new  = amrex::ReduceSum(pc_new, f);
        ParallelDescriptor::ReduceLongSum(sm_orig);
        ParallelDescriptor::ReduceLongSum(sm_new);
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(sm_orig == sm_new,
            "Int component sum mismatch after restart (comp " + std::to_string(icomp) + ")");
    }
}

void test ()
{
    const int ncells        = 32;
    const int max_grid_size = 16;
    const int nppc          = 2;
    const int iseed         = 451;

    Box domain(IntVect(0), IntVect(ncells-1));
    BoxArray ba(domain);
    ba.maxSize(max_grid_size);

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    Array<int,AMREX_SPACEDIM> is_per{AMREX_D_DECL(1,1,1)};
    Geometry geom(domain, real_box, CoordSys::cartesian, is_per);
    DistributionMapping dm(ba);

    MyPC pc(geom, dm, ba);
    pc.SetVerbose(false);
    MyPC::ParticleInitData pdata = {{1.0, 2.0}, {3}, {4.0, 5.0}, {6}};
    pc.InitRandom(nppc * AMREX_D_TERM(ncells, *ncells, *ncells),
                  iseed, pdata, /*serialize=*/false);

    Vector<std::string> real_names, int_names;
    for (int i = 0; i < NStructReal + NArrayReal; ++i) {
        real_names.push_back("real_" + std::to_string(i));
    }
    for (int i = 0; i < NStructInt + NArrayInt; ++i) {
        int_names.push_back("int_" + std::to_string(i));
    }

    // Reference checkpoint without pre/post.
    pc.CheckpointHDF5("chk_ref", "particles", true, real_names, int_names);

    // Same data written through the pre/post path.
    pc.SetUsePrePost(true);
    pc.CheckpointPreHDF5();
    pc.CheckpointHDF5("chk_prepost", "particles", true, real_names, int_names);
    pc.CheckpointPostHDF5();
    pc.SetUsePrePost(false);

    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor()) {
        auto ref_size     = std::filesystem::file_size("chk_ref/particles/particles.h5");
        auto prepost_size = std::filesystem::file_size("chk_prepost/particles/particles.h5");
        amrex::Print() << "  reference HDF5 file size: " << ref_size
                       << ", pre/post HDF5 file size: " << prepost_size << "\n";
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(ref_size == prepost_size,
            "CheckpointPostHDF5 modified the HDF5 file");
    }

    MyPC pc_read(geom, dm, ba);
    pc_read.SetVerbose(false);
    pc_read.RestartHDF5("chk_prepost/particles", "particles");

    verify_same(pc, pc_read);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    test();

    amrex::Print() << "HDF5 pre/post checkpoint/restart test PASSED\n";

    amrex::Finalize();
}
