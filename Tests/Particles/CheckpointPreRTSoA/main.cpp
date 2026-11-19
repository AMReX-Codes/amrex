#include <AMReX.H>
#include <AMReX_ParticleContainer.H>

using namespace amrex;
using PC = ParticleContainerRTSoA<>;

void schema (PC& pc)
{
    pc.SetArena(The_Arena());
    for (int d = 0; d < AMREX_SPACEDIM; ++d) { pc.AddRealComp("pos" + std::to_string(d)); }
    pc.AddRealComp("mass");
    pc.AddRealComp("tag");
    pc.AddIntComp("first");
    pc.AddIntComp("second");
}

int main (int argc, char* argv[])
{
    Initialize(argc, argv);
    {
        amrex::RealBox bounds;
        for (int d = 0; d < AMREX_SPACEDIM; ++d) { bounds.setLo(d,0); bounds.setHi(d,1); }
        int periodic[AMREX_SPACEDIM] = {AMREX_D_DECL(0,0,0)};
        Vector<Geometry> geom;
        geom.emplace_back(Box(IntVect(0), IntVect(7)), &bounds, 0, periodic);
        geom.emplace_back(Box(IntVect(0), IntVect(15)), &bounds, 0, periodic);
        Vector<BoxArray> ba{BoxArray(geom[0].Domain()), BoxArray(Box(IntVect(0), IntVect(7)))};
        for (auto& boxes : ba) { boxes.maxSize(4); }
        Vector<DistributionMapping> dm{DistributionMapping(ba[0]), DistributionMapping(ba[1])};
        Vector<IntVect> rr{IntVect(2)};
        Vector<Long> expected(2, 0);
        for (int lev = 0; lev < 2; ++lev) {
            for (int grid = 0; grid < ba[lev].size(); ++grid) {
                if (lev == 1 || ba[lev][grid].smallEnd(0) >= 4) { ++expected[lev]; }
            }
        }
        for (bool prepost : {false, true}) {
            PC src(geom, dm, ba, rr);
            schema(src);
            src.SetUsePrePost(prepost); // Isolate this test from instance-option parsing.
            for (int lev = 0; lev < 2; ++lev) {
                for (MFIter mfi = src.MakeMFIter(lev); mfi.isValid(); ++mfi) {
                    auto& tile = src.DefineAndReturnParticleTile(lev, mfi);
                    if (lev == 0 && mfi.validbox().smallEnd(0) < 4) { continue; }
                    tile.resize(2);
                    auto ptd = tile.getParticleTileData();
                    auto dx = geom[lev].CellSizeArray();
                    auto box = mfi.validbox();
                    Long id = 10000*(lev+1) + mfi.index() + 1;
                    ParallelFor(2, [=] AMREX_GPU_DEVICE (int i) noexcept {
                        ptd.id(i) = (i == 0) ? id : -id;
                        ptd.cpu(i) = 0;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            ptd.pos(d,i) = (Real(box.smallEnd(d)) + Real(0.5))*dx[d];
                        }
                        ptd.rdata(AMREX_SPACEDIM,i) = ParticleReal(id);
                        ptd.rdata(AMREX_SPACEDIM+1,i) = ParticleReal(id)*2;
                        ptd.idata(0,i) = int(id)+7;
                        ptd.idata(1,i) = -int(id);
                    });
                }
            }
            Gpu::streamSynchronize();
            PC::ParticleType::NextID(40000 + ParallelDescriptor::MyProc());
            auto dir = prepost ? "checkpoint_pre" : "checkpoint_normal";
            src.CheckpointPre();
            Vector<std::string> real_names{"mass", "tag"};
            Vector<std::string> int_names{"first", "second"};
            src.Checkpoint(dir, "particles", true, real_names, int_names);
            src.CheckpointPost();
            auto header = ParticleHeader::read(dir, "particles");
            AMREX_ALWAYS_ASSERT(header.num_particles == expected[0]+expected[1]);
            AMREX_ALWAYS_ASSERT(header.next_id == 40000 + ParallelDescriptor::NProcs()-1);
            AMREX_ALWAYS_ASSERT(header.real_comp_names == real_names);
            AMREX_ALWAYS_ASSERT(header.int_comp_names == int_names);
            AMREX_ALWAYS_ASSERT(header.finest_level == 1);
            for (int lev = 0; lev < 2; ++lev) {
                Long n = 0;
                for (auto const& entry : header.grids[lev]) { n += entry.count; }
                AMREX_ALWAYS_ASSERT(n == expected[lev]);
            }
            PC dst(geom, dm, ba, rr);
            schema(dst);
            dst.Restart(dir, "particles");
            Long next = PC::ParticleType::NextID();
            AMREX_ALWAYS_ASSERT(next == 40000 + ParallelDescriptor::NProcs()-1);
            PC::ParticleType::NextID(next);
            for (int lev = 0; lev < 2; ++lev) {
                AMREX_ALWAYS_ASSERT(dst.NumberOfParticlesAtLevel(lev) == expected[lev]);
                ReduceOps<ReduceOpSum> op;
                ReduceData<int> data(op);
                using Tuple = decltype(data)::Type;
                for (auto const& kv : dst.GetParticles(lev)) {
                    auto ptd = kv.second.getConstParticleTileData();
                    int grid = kv.first.first;
                    auto box = ba[lev][grid];
                    auto dx = geom[lev].CellSizeArray();
                    Long expected_id = 10000*(lev+1) + grid + 1;
                    op.eval(kv.second.numParticles(), data, [=] AMREX_GPU_DEVICE (int i) -> Tuple {
                        int bad = ptd.id(i) != expected_id || ptd.cpu(i) != 0;
                        for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                            bad += ptd.pos(d,i) != ParticleReal((Real(box.smallEnd(d))+Real(0.5))*dx[d]);
                        }
                        bad += ptd.rdata(AMREX_SPACEDIM,i) != ParticleReal(expected_id);
                        bad += ptd.rdata(AMREX_SPACEDIM+1,i) != ParticleReal(2*expected_id);
                        bad += ptd.idata(0,i) != expected_id+7;
                        bad += ptd.idata(1,i) != -expected_id;
                        return {bad};
                    });
                }
                AMREX_ALWAYS_ASSERT(amrex::get<0>(data.value(op)) == 0);
            }
        }
    }
    Finalize();
}
