#include <AMReX.H>
#include <AMReX_ParticleContainer.H>
#include <AMReX_Print.H>

#include <fstream>

using namespace amrex;
using PC = ParticleContainerRTSoA<>;

void check (PC const& pc, Long count, Long first_id, bool ascii, bool defaults, bool cells)
{
    AMREX_ALWAYS_ASSERT(pc.TotalNumberOfParticles() == count);
    Long next = PC::ParticleType::NextID();
    PC::ParticleType::NextID(next);
    Long allocated = next-first_id;
    ParallelDescriptor::ReduceLongSum(allocated);
    AMREX_ALWAYS_ASSERT(allocated == count);
    Vector<int> ids(count*ParallelDescriptor::NProcs(),0);
    Vector<int> cell_count(pc.Geom(0).Domain().numPts(),0);
    for (auto const& kv : pc.GetParticles(0)) {
        auto const& tile = kv.second;
        AMREX_ALWAYS_ASSERT(tile.arena() == pc.arena());
        AMREX_ALWAYS_ASSERT(tile.GetRealData("mass").data() == tile.GetRealData(AMREX_SPACEDIM).data());
        AMREX_ALWAYS_ASSERT(tile.GetIntData("first").data() == tile.GetIntData(0).data());
        auto n = tile.numParticles();
        Gpu::HostVector<uint64_t> identity(n);
        Gpu::copy(Gpu::deviceToHost, tile.GetIdCPUData().begin(), tile.GetIdCPUData().begin()+n, identity.begin());
        Vector<Gpu::HostVector<ParticleReal>> real(pc.NumRealComps());
        Vector<Gpu::HostVector<int>> integer(pc.NumIntComps());
        for (int c = 0; c < pc.NumRealComps(); ++c) {
            real[c].resize(n);
            Gpu::copy(Gpu::deviceToHost, tile.GetRealData(c).begin(), tile.GetRealData(c).begin()+n, real[c].begin());
        }
        for (int c = 0; c < pc.NumIntComps(); ++c) {
            integer[c].resize(n);
            Gpu::copy(Gpu::deviceToHost, tile.GetIntData(c).begin(), tile.GetIntData(c).begin()+n, integer[c].begin());
        }
        for (Long i = 0; i < n; ++i) {
            Long id = ParticleIDWrapper(identity[i]);
            int cpu = ParticleCPUWrapper(identity[i]);
            AMREX_ALWAYS_ASSERT(id >= first_id && id < first_id+count);
            AMREX_ALWAYS_ASSERT(cpu >= 0 && cpu < ParallelDescriptor::NProcs());
            ++ids[cpu*count+id-first_id];
            Long cell_index = 0;
            Long stride = 1;
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                auto x = real[d][i];
                AMREX_ALWAYS_ASSERT(x >= 0 && x < 1);
                if (cells) {
                    auto scaled = x*8;
                    int cell = int(scaled);
                    AMREX_ALWAYS_ASSERT(scaled == ParticleReal(cell)+ParticleReal(0.5));
                    cell_index += cell*stride;
                    stride *= 8;
                }
            }
            if (cells) { ++cell_count[cell_index]; }
            auto mass = real[AMREX_SPACEDIM][i];
            AMREX_ALWAYS_ASSERT(ascii ? (mass == 7 || mass == 9) : mass == (defaults ? 0 : 3.25));
            AMREX_ALWAYS_ASSERT(real[AMREX_SPACEDIM+1][i] == (ascii || defaults ? 0 : -6.5));
            AMREX_ALWAYS_ASSERT(integer[0][i] == (ascii || defaults ? 0 : 17));
            AMREX_ALWAYS_ASSERT(integer[1][i] == (ascii || defaults ? 0 : -23));
        }
    }
    ParallelDescriptor::ReduceIntSum(ids.data(), ids.size());
    for (int n : ids) { AMREX_ALWAYS_ASSERT(n <= 1); }
    if (cells) {
        ParallelDescriptor::ReduceIntSum(cell_count.data(), cell_count.size());
        for (int n : cell_count) { AMREX_ALWAYS_ASSERT(n == 1); }
    }
}

int main (int argc, char* argv[])
{
    Initialize(argc, argv);
    {
        // A nonzero integer domain origin checks the index-to-position conversion.
        Box domain(IntVect(4), IntVect(11));
        RealBox bounds;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { bounds.setLo(d,0); bounds.setHi(d,1); }
        int periodic[AMREX_SPACEDIM] = {AMREX_D_DECL(0,0,0)};
        Geometry geom(domain,&bounds,0,periodic);
        BoxArray ba(domain);
        ba.maxSize(4);
        DistributionMapping dm(ba);
        if (ParallelDescriptor::IOProcessor()) {
            std::ofstream file("particles.txt");
            file << "2\n";
            for (int p = 0; p < 2; ++p) {
                for (int d = 0; d < AMREX_SPACEDIM; ++d) { file << (p == 0 ? 0.125 : 0.25) << ' '; }
                file << (p == 0 ? 7 : 9);
                if (p == 0) { file << '\n'; }
            }
        }
        ParallelDescriptor::Barrier();
        for (int mode = 0; mode < 5; ++mode) {
            PC pc(geom,dm,ba);
            pc.SetArena(The_Arena());
            for (int d = 0; d < AMREX_SPACEDIM; ++d) { pc.AddRealComp("pos"+std::to_string(d)); }
            pc.AddRealComp("mass"); pc.AddRealComp("tag");
            pc.AddIntComp("first"); pc.AddIntComp("second");
            // Predefine empty tiles with poisoned capacity to expose missing writes.
            for (MFIter mfi = pc.MakeMFIter(0); mfi.isValid(); ++mfi) {
                if (mfi.index()%2 != 0) { continue; }
                auto& tile = pc.DefineAndReturnParticleTile(0,mfi);
                tile.resize(32);
                auto ptd = tile.getParticleTileData();
                ParallelFor(32, [=] AMREX_GPU_DEVICE (int i) noexcept {
                    for (int c = 0; c < AMREX_SPACEDIM+2; ++c) { ptd.rdata(c,i) = -999; }
                    for (int c = 0; c < 2; ++c) { ptd.idata(c,i) = -999; }
                });
                Gpu::streamSynchronize();
                tile.resize(0);
            }
            PC::ParticleType::NextID(1000);
            PC::ParticleInitData data{{3.25,-6.5},{17,-23}};
            Long count = 37; // Uneven MPI work partition.
            if (mode == 0) {
                IntVect reps(2);
                pc.InitFromAsciiFile("particles.txt",1,&reps);
                count = 2*(1 << AMREX_SPACEDIM);
            } else if (mode == 3) {
                pc.InitOnePerCell(0.5,0.5,0.5,data);
                count = domain.numPts();
            } else {
                if (mode == 4) { data = {}; }
                pc.InitRandom(count,12345,data,mode == 1);
                if (mode == 1 || mode == 2) {
                    pc.InitRandom(13,54321,data,mode == 1);
                    count += 13;
                }
            }
            check(pc,count,1000,mode == 0,mode == 4,mode == 3);
            pc.Redistribute();
            check(pc,count,1000,mode == 0,mode == 4,mode == 3);
        }
        Print() << "RTSoA initializer counts, values, defaults, positions and IDs passed\n";
    }
    Finalize();
}
