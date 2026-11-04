// Redistribute() with a custom cell assignor that does not take the last direction's cell
// from the position, as for particles on a terrain-following grid that carry their vertical
// cell index while their position holds the physical height. Such a particle can lie outside
// the problem domain in that direction and still belong to a cell inside it; it must be kept,
// as on the GPU path, on every Redistribute path: the full one with several boxes, the
// in-place one with a single box, and across a periodic boundary in the other directions.
// A particle whose stored cell is outside the domain must still be removed.
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Particles.H>
#include <AMReX_Print.H>

#include <string>

using namespace amrex;

namespace
{

constexpr int last_dir = AMREX_SPACEDIM - 1;

//! The cell from the position, except in the last direction, where it is idata(0)
struct StoredCellAssignor
{
    template <typename P>
    AMREX_GPU_HOST_DEVICE
    IntVect operator() (P const& p,
                        GpuArray<Real,AMREX_SPACEDIM> const& plo,
                        GpuArray<Real,AMREX_SPACEDIM> const& dxi,
                        const Box& domain) const noexcept
    {
        IntVect iv = getParticleCell(p, plo, dxi, domain);
        iv[last_dir] = p.idata(0) + domain.smallEnd(last_dir);
        return iv;
    }
};

using PC = ParticleContainer<0, 1, 0, 0, DefaultAllocator, StoredCellAssignor>;

enum struct Layout { OneBox, Boxes };

struct Counts
{
    Long n;          //!< valid particles
    Long n_kept;     //!< of those, the ones whose stored cell is inside the domain
};

void test (Layout layout, bool periodic)
{
    int const n_cell = 16;
    Box const domain(IntVect(0), IntVect(n_cell - 1));
    RealBox const real_box({AMREX_D_DECL(Real(0.), Real(0.), Real(0.))},
                           {AMREX_D_DECL(Real(1.), Real(1.), Real(1.))});
    // periodic, if at all, in every direction but the last
    Array<int, AMREX_SPACEDIM> is_per{AMREX_D_DECL(int(periodic), int(periodic), 0)};
    is_per[last_dir] = 0;
    Geometry const geom(domain, real_box, CoordSys::cartesian, is_per);

    BoxArray ba(domain);
    if (layout == Layout::Boxes) { ba.maxSize(n_cell / 4); }
    DistributionMapping const dm(ba);

    PC pc(geom, dm, ba);

    // one particle at the center of every cell, with its vertical cell stored
    {
        using HostTile = ParticleTile<PC::ParticleType, 0, 0, PinnedArenaAllocator>;
        auto const dx = geom.CellSizeArray();
        auto const plo = geom.ProbLoArray();
        for (MFIter mfi = pc.MakeMFIter(0); mfi.isValid(); ++mfi) {
            HostTile host_tile;
            host_tile.define(0, 0);
            Box const& bx = mfi.tilebox();
            for (IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
                PC::ParticleType p;
                p.id() = PC::ParticleType::NextID();
                p.cpu() = ParallelDescriptor::MyProc();
                for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                    p.pos(d) = ParticleReal(plo[d] + (Real(iv[d]) + Real(0.5)) * dx[d]);
                }
                p.idata(0) = iv[last_dir];
                host_tile.push_back(p);
            }
            auto& tile = pc.DefineAndReturnParticleTile(0, mfi);
            tile.resize(host_tile.numParticles());
            amrex::copyParticles(tile, host_tile);
        }
        Gpu::streamSynchronize();
    }
    AMREX_ALWAYS_ASSERT(pc.OK());

    // Move them: the last direction's position out of the domain for all (the stored cell
    // still places them), some across the upper boundary of the first direction (periodic or
    // not), some into the next box, and a few with a stored cell outside the domain.
    auto const dx = geom.CellSizeArray();
    auto const rhi = geom.ProbHiArrayInParticleReal();
    auto const rlo = geom.ProbLoArrayInParticleReal();
    ReduceOps<ReduceOpSum, ReduceOpSum> reduce_op;
    ReduceData<Long, Long> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    for (auto& kv : pc.GetParticles(0)) {
        auto& ptile = kv.second;
        auto* pstruct = ptile.GetArrayOfStructs()().data();
        reduce_op.eval(ptile.numParticles(), reduce_data,
        [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
        {
            auto& p = pstruct[i];
            Long const pid = p.id();
            p.pos(last_dir) = (pid % 2 == 0) ? rhi[last_dir] + ParticleReal(0.25)
                                             : rlo[last_dir] - ParticleReal(0.25);
            if (pid % 5 == 0) {
                p.pos(0) += rhi[0] - rlo[0];                        // one domain length on
            } else if (pid % 5 == 1) {
                p.pos(0) = amrex::min(p.pos(0) + ParticleReal(4 * dx[0]),
                                      rhi[0] - ParticleReal(0.5 * dx[0]));  // next box
            }
            bool kept = true;
            if (pid % 7 == 0) {
                p.idata(0) = n_cell;                                 // stored cell outside
                kept = false;
            }
            if (pid % 5 == 0 && !periodic) { kept = false; }         // left the domain
            return {1, kept ? 1 : 0};
        });
    }
    auto const hv = reduce_data.value(reduce_op);
    Counts c{amrex::get<0>(hv), amrex::get<1>(hv)};
    ParallelDescriptor::ReduceLongSum(c.n);
    ParallelDescriptor::ReduceLongSum(c.n_kept);

    pc.Redistribute();

    Long const n_after = pc.TotalNumberOfParticles();
    std::string const name = std::string(layout == Layout::OneBox ? "one box" : "boxes")
                             + (periodic ? ", periodic" : ", not periodic");
    amrex::Print() << name << ": " << c.n << " particles, " << c.n_kept
                   << " expected to stay, " << n_after << " after Redistribute\n";
    AMREX_ALWAYS_ASSERT(n_after == c.n_kept);
    AMREX_ALWAYS_ASSERT(pc.OK());

    // the kept particles are back inside in the first direction and untouched in the last
    ReduceOps<ReduceOpSum> check_op;
    ReduceData<Long> check_data(check_op);
    for (auto& kv : pc.GetParticles(0)) {
        auto const* pstruct = kv.second.GetArrayOfStructs()().data();
        check_op.eval(kv.second.numParticles(), check_data,
        [=] AMREX_GPU_DEVICE (int i) -> GpuTuple<Long>
        {
            auto const& p = pstruct[i];
            bool const ok = p.pos(0) >= rlo[0] && p.pos(0) <= rhi[0]
                && (p.pos(last_dir) == rhi[last_dir] + ParticleReal(0.25)
                    || p.pos(last_dir) == rlo[last_dir] - ParticleReal(0.25));
            return {ok ? 0 : 1};
        });
    }
    Long n_bad = amrex::get<0>(check_data.value(check_op));
    ParallelDescriptor::ReduceLongSum(n_bad);
    AMREX_ALWAYS_ASSERT(n_bad == 0);
}

}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        for (bool periodic : {true, false}) {
            test(Layout::OneBox, periodic);
            test(Layout::Boxes, periodic);
        }
        amrex::Print() << "PASSED\n";
    }
    amrex::Finalize();
}
