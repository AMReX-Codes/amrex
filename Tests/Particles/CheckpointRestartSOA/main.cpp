/**
 * Test for dual-grid particle checkpoint/restart capability using pure SoA
 * particles and HDF5 I/O.
 *
 * Verifies that particles are preserved correctly when reading a
 * checkpoint file into a ParticleContainerPureSoA with a different mesh structure,
 * covering:
 *   1. Restart with fewer AMR levels than the checkpoint
 *   2. Restart with more AMR levels than the checkpoint
 *   3. Restart with the same number of levels but different BoxArrays
 */
#include <AMReX.H>
#include <AMReX_Particles.H>

#include <algorithm>
#include <array>
#include <cstring>

using namespace amrex;

// Pure SoA particle type: 6 real components (including positions), 2 int components
constexpr int NReal = 6;
constexpr int NInt  = 2;

using MyPC        = ParticleContainerPureSoA<NReal, NInt>;
struct MeshData {
    Vector<Box>                 domains;
    Vector<BoxArray>            ba;
    Vector<IntVect>             ref_ratio;
    Vector<Geometry>            geom;
    Vector<DistributionMapping> dmap;
};

struct ParticleRecord {
    std::uint64_t idcpu = 0;
    std::array<ParticleReal, NReal> real{};
    std::array<int, NInt> idata{};

    bool operator< (ParticleRecord const& rhs) const
    {
        if (idcpu != rhs.idcpu) { return idcpu < rhs.idcpu; }
        if (real != rhs.real) { return real < rhs.real; }
        return idata < rhs.idata;
    }

    bool operator== (ParticleRecord const& rhs) const
    {
        return idcpu == rhs.idcpu && real == rhs.real && idata == rhs.idata;
    }
};

Vector<ParticleRecord>
collect_local_records (MyPC const& pc)
{
    Vector<ParticleRecord> records;

    for (int lev = 0; lev <= pc.finestLevel(); ++lev) {
        for (auto const& kv : pc.GetParticles(lev)) {
            auto const& ptile = kv.second;
            auto const np = ptile.numParticles();
            auto const ptd = ptile.getConstParticleTileData();

            for (int i = 0; i < np; ++i) {
                ParticleRecord rec;
                rec.idcpu = ptd.idcpu(i);
                for (int icomp = 0; icomp < NReal; ++icomp) {
                    rec.real[icomp] = ptd.rdata(icomp)[i];
                }
                for (int icomp = 0; icomp < NInt; ++icomp) {
                    rec.idata[icomp] = ptd.idata(icomp)[i];
                }
                records.push_back(rec);
            }
        }
    }

    return records;
}

Vector<ParticleRecord>
gather_sorted_records (MyPC const& pc)
{
    auto local_records = collect_local_records(pc);
    auto const local_nbytes = static_cast<int>(local_records.size() * sizeof(ParticleRecord));
    auto const ioproc = ParallelDescriptor::IOProcessorNumber();
    auto const recv_counts = ParallelDescriptor::Gather(local_nbytes, ioproc);

    Vector<int> recv_offsets;
    Vector<char> recv_buffer;

    if (ParallelDescriptor::IOProcessor()) {
        recv_offsets.resize(recv_counts.size(), 0);
        int total_nbytes = 0;
        for (int i = 0, n = static_cast<int>(recv_counts.size()); i < n; ++i) {
            recv_offsets[i] = total_nbytes;
            total_nbytes += recv_counts[i];
        }
        recv_buffer.resize(total_nbytes);
    }

    auto const* send_ptr = local_records.empty()
        ? nullptr
        : reinterpret_cast<char const*>(local_records.data());
    auto* recv_ptr = recv_buffer.empty() ? nullptr : recv_buffer.data();

    ParallelDescriptor::Gatherv(send_ptr, local_nbytes, recv_ptr,
                                recv_counts, recv_offsets, ioproc);

    Vector<ParticleRecord> gathered_records;
    if (ParallelDescriptor::IOProcessor()) {
        gathered_records.resize(recv_buffer.size() / sizeof(ParticleRecord));
        if (!recv_buffer.empty()) {
            std::memcpy(gathered_records.data(), recv_buffer.data(), recv_buffer.size());
        }
        std::sort(gathered_records.begin(), gathered_records.end());
    }

    return gathered_records;
}

/**
 * Build a nested AMR hierarchy with `nlevs` levels over a [0,1]^d domain
 * discretised by `ncells` cells in each direction. The coarse BoxArray is
 * decomposed into boxes of at most `max_grid_size` cells; the fine level
 * covers the central half of the domain.
 */
MeshData build_mesh (int ncells, int nlevs, int max_grid_size)
{
    AMREX_ALWAYS_ASSERT(nlevs >= 1 && nlevs <= 2);

    MeshData m;

    // Level-0 domain
    IntVect lo(AMREX_D_DECL(0, 0, 0));
    IntVect hi(AMREX_D_DECL(ncells-1, ncells-1, ncells-1));

    m.domains.resize(nlevs);
    m.domains[0].setSmall(lo);
    m.domains[0].setBig(hi);

    m.ref_ratio.resize(nlevs > 1 ? nlevs - 1 : 0);
    for (int lev = 1; lev < nlevs; ++lev) {
        m.ref_ratio[lev-1] = IntVect(AMREX_D_DECL(2, 2, 2));
        m.domains[lev] = amrex::refine(m.domains[lev-1], m.ref_ratio[lev-1]);
    }

    m.ba.resize(nlevs);
    m.ba[0].define(m.domains[0]);
    m.ba[0].maxSize(max_grid_size);

    if (nlevs > 1) {
        // Refined region: the central 1/4 of the fine-level domain
        int n_fine = ncells * 2; // ref_ratio == 2
        IntVect rlo(AMREX_D_DECL(n_fine/4,   n_fine/4,   n_fine/4));
        IntVect rhi(AMREX_D_DECL(3*n_fine/4-1, 3*n_fine/4-1, 3*n_fine/4-1));
        m.ba[1].define(Box(rlo, rhi));
        m.ba[1].maxSize(max_grid_size);
    }

    RealBox real_box;
    for (int n = 0; n < AMREX_SPACEDIM; ++n) {
        real_box.setLo(n, 0.0);
        real_box.setHi(n, 1.0);
    }
    int is_per[] = {AMREX_D_DECL(1, 1, 1)};

    m.geom.resize(nlevs);
    m.geom[0].define(m.domains[0], &real_box, CoordSys::cartesian, is_per);
    for (int lev = 1; lev < nlevs; ++lev) {
        m.geom[lev].define(m.domains[lev], &real_box, CoordSys::cartesian, is_per);
    }

    m.dmap.resize(nlevs);
    for (int lev = 0; lev < nlevs; ++lev) {
        m.dmap[lev] = DistributionMapping{m.ba[lev]};
    }

    return m;
}

/**
 * Assert that two ParticleContainerPureSoA hold identical particle data by
 * comparing exact sorted particle records gathered on the I/O rank.
 */
void verify_same (MyPC& pc_orig, MyPC& pc_new)
{
    auto orig_records = gather_sorted_records(pc_orig);
    auto new_records = gather_sorted_records(pc_new);

    if (ParallelDescriptor::IOProcessor()) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            orig_records.size() == new_records.size(),
            "Particle count mismatch after restart");

        for (Long i = 0; i < orig_records.size(); ++i) {
            AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
                orig_records[i] == new_records[i],
                "Particle data mismatch after restart");
        }
    }
}

/**
 * Test 0: Write and restart with the same hierarchy and BoxArray.
 * This preserves the original checkpoint/restart coverage without any
 * hierarchy or decomposition changes.
 */
void test_same_hierarchy ()
{
    amrex::Print() << "Test 0: restart with same hierarchy\n";

    const int ncells         = 32;
    const int max_grid_size  = 16;
    const int nppc           = 2;
    const int iseed          = 451;
    const std::string chkdir = "chk_same_hierarchy";

    MyPC::ParticleInitData pdata = {{}, {}, {1.0, 2.0, 3.0, 4.0, 6.0, 7.0}, {5, 8}};

    Vector<std::string> real_names, int_names;
    for (int i = 0; i < NReal - AMREX_SPACEDIM; ++i) {
        real_names.push_back("real_" + std::to_string(i));
    }
    for (int i = 0; i < NInt; ++i) {
        int_names.push_back("int_" + std::to_string(i));
    }

    auto mesh = build_mesh(ncells, 2, max_grid_size);
    MyPC pc_write(mesh.geom, mesh.dmap, mesh.ba, mesh.ref_ratio);
    pc_write.SetVerbose(false);
    pc_write.InitRandom(nppc * AMREX_D_TERM(ncells, *ncells, *ncells),
                        iseed, pdata, /*serialize=*/false);
#ifdef AMREX_USE_HDF5
    pc_write.CheckpointHDF5(chkdir, "particles", true, real_names, int_names);
#else
    pc_write.Checkpoint(chkdir, "particles", real_names, int_names);
#endif

    AsyncOut::Finish();
    ParallelDescriptor::Barrier();

    MyPC pc_read(mesh.geom, mesh.dmap, mesh.ba, mesh.ref_ratio);
    pc_read.SetVerbose(false);
#ifdef AMREX_USE_HDF5
    pc_read.RestartHDF5(chkdir + "/particles", "particles");
#else
    pc_read.Restart(chkdir, "particles");
#endif

    verify_same(pc_write, pc_read);

    amrex::Print() << "  PASSED\n";
}

/**
 * Test 1: Write a 2-level HDF5 checkpoint, restart into a 1-level container.
 * Particles that resided on the fine level must be redistributed to the
 * coarse level; totals and component sums must be unchanged.
 */
void test_fewer_levels ()
{
    amrex::Print() << "Test 1: restart with fewer levels than checkpoint\n";

    const int ncells         = 32;
    const int max_grid_size  = 16;
    const int nppc           = 2;
    const int iseed          = 451;
    const std::string chkdir = "chk_fewer_levels";

    MyPC::ParticleInitData pdata = {{}, {}, {1.0, 2.0, 3.0, 4.0, 6.0, 7.0}, {5, 8}};

    Vector<std::string> real_names, int_names;
    for (int i = 0; i < NReal - AMREX_SPACEDIM; ++i) {
        real_names.push_back("real_" + std::to_string(i));
    }
    for (int i = 0; i < NInt; ++i) {
        int_names.push_back("int_" + std::to_string(i));
    }

    // --- write with 2 levels ---
    auto mesh2 = build_mesh(ncells, 2, max_grid_size);
    MyPC pc_write(mesh2.geom, mesh2.dmap, mesh2.ba, mesh2.ref_ratio);
    pc_write.SetVerbose(false);
    pc_write.InitRandom(nppc * AMREX_D_TERM(ncells, *ncells, *ncells),
                        iseed, pdata, /*serialize=*/false);
#ifdef AMREX_USE_HDF5
    pc_write.CheckpointHDF5(chkdir, "particles", true, real_names, int_names);
#else
    pc_write.Checkpoint(chkdir, "particles", real_names, int_names);
#endif

    AsyncOut::Finish();
    ParallelDescriptor::Barrier();

    // --- restart with 1 level ---
    auto mesh1 = build_mesh(ncells, 1, max_grid_size);
    MyPC pc_read(mesh1.geom, mesh1.dmap, mesh1.ba, mesh1.ref_ratio);
    pc_read.SetVerbose(false);
#ifdef AMREX_USE_HDF5
    pc_read.RestartHDF5(chkdir + "/particles", "particles");
#else
    pc_read.Restart(chkdir, "particles");
#endif

    verify_same(pc_write, pc_read);

    amrex::Print() << "  PASSED\n";
}

/**
 * Test 2: Write a 1-level HDF5 checkpoint, restart into a 2-level container.
 * After restart, Redistribute() assigns particles inside the refined
 * region to the finer level; totals and component sums must be unchanged.
 */
void test_more_levels ()
{
    amrex::Print() << "Test 2: restart with more levels than checkpoint\n";

    const int ncells         = 32;
    const int max_grid_size  = 16;
    const int nppc           = 2;
    const int iseed          = 451;
    const std::string chkdir = "chk_more_levels";

    MyPC::ParticleInitData pdata = {{}, {}, {1.0, 2.0, 3.0, 4.0, 6.0, 7.0}, {5, 8}};

    Vector<std::string> real_names, int_names;
    for (int i = 0; i < NReal - AMREX_SPACEDIM; ++i) {
        real_names.push_back("real_" + std::to_string(i));
    }
    for (int i = 0; i < NInt; ++i) {
        int_names.push_back("int_" + std::to_string(i));
    }

    // --- write with 1 level ---
    auto mesh1 = build_mesh(ncells, 1, max_grid_size);
    MyPC pc_write(mesh1.geom, mesh1.dmap, mesh1.ba, mesh1.ref_ratio);
    pc_write.SetVerbose(false);
    pc_write.InitRandom(nppc * AMREX_D_TERM(ncells, *ncells, *ncells),
                        iseed, pdata, /*serialize=*/false);
#ifdef AMREX_USE_HDF5
    pc_write.CheckpointHDF5(chkdir, "particles", true, real_names, int_names);
#else
    pc_write.Checkpoint(chkdir, "particles", real_names, int_names);
#endif

    AsyncOut::Finish();
    ParallelDescriptor::Barrier();

    // --- restart with 2 levels ---
    auto mesh2 = build_mesh(ncells, 2, max_grid_size);
    MyPC pc_read(mesh2.geom, mesh2.dmap, mesh2.ba, mesh2.ref_ratio);
    pc_read.SetVerbose(false);
#ifdef AMREX_USE_HDF5
    pc_read.RestartHDF5(chkdir + "/particles", "particles");
#else
    pc_read.Restart(chkdir, "particles");
#endif

    verify_same(pc_write, pc_read);

    amrex::Print() << "  PASSED\n";
}

/**
 * Test 3: Write an HDF5 checkpoint with one BoxArray decomposition, restart
 * into a container with the same number of levels but a different BoxArray
 * (different max_grid_size). This is the canonical dual-grid scenario: the
 * Particle_H header written at checkpoint time differs from the BoxArray
 * of the new container, so AMReX uses temporary grids while reading and
 * calls Redistribute() to move particles onto the new decomposition.
 */
void test_different_boxarrays ()
{
    amrex::Print() << "Test 3: restart with same levels but different BoxArrays\n";

    const int ncells         = 32;
    const int nppc           = 2;
    const int iseed          = 451;
    const std::string chkdir = "chk_different_ba";

    MyPC::ParticleInitData pdata = {{}, {}, {1.0, 2.0, 3.0, 4.0, 6.0, 7.0}, {5, 8}};

    Vector<std::string> real_names, int_names;
    for (int i = 0; i < NReal - AMREX_SPACEDIM; ++i) {
        real_names.push_back("real_" + std::to_string(i));
    }
    for (int i = 0; i < NInt; ++i) {
        int_names.push_back("int_" + std::to_string(i));
    }

    // --- write: coarse decomposition (few large boxes) ---
    auto mesh_coarse = build_mesh(ncells, 1, ncells); // one box per rank at most
    MyPC pc_write(mesh_coarse.geom, mesh_coarse.dmap, mesh_coarse.ba,
                  mesh_coarse.ref_ratio);
    pc_write.SetVerbose(false);
    pc_write.InitRandom(nppc * AMREX_D_TERM(ncells, *ncells, *ncells),
                        iseed, pdata, /*serialize=*/false);
#ifdef AMREX_USE_HDF5
    pc_write.CheckpointHDF5(chkdir, "particles", true, real_names, int_names);
#else
    pc_write.Checkpoint(chkdir, "particles", real_names, int_names);
#endif

    AsyncOut::Finish();
    ParallelDescriptor::Barrier();

    // --- restart: finer decomposition (many smaller boxes) ---
    auto mesh_fine = build_mesh(ncells, 1, ncells / 4);
    MyPC pc_read(mesh_fine.geom, mesh_fine.dmap, mesh_fine.ba,
                 mesh_fine.ref_ratio);
    pc_read.SetVerbose(false);
#ifdef AMREX_USE_HDF5
    pc_read.RestartHDF5(chkdir + "/particles", "particles");
#else
    pc_read.Restart(chkdir, "particles");
#endif

    verify_same(pc_write, pc_read);

    amrex::Print() << "  PASSED\n";
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);

    test_same_hierarchy();
    test_fewer_levels();
    test_more_levels();
    test_different_boxarrays();

    amrex::Print() << "All dual-grid HDF5 SOA restart tests PASSED\n";

    amrex::Finalize();
}
