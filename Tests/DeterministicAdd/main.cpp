// Test deterministic add mode by comparing with standard SumBoundary
#include <AMReX.H>
#include <AMReX_BoxArray.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_MultiFab.H>
#include <AMReX_FabArray.H>
#include <AMReX_Geometry.H>
#include <random>

using namespace amrex;

static void fill_values (MultiFab& mf)
{
    // Set different constants on each local fab to create distinct contributions
    int idx = 0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi, ++idx) {
        Real val = Real(1 + idx); // 1, 2, ...
        mf[mfi].setVal(val);
    }
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        int iseed = 1; amrex::ignore_unused(iseed);
        int n_cell = 16;
        int max_grid_size = 8;
        int ncomp = 1;
        int ng = 1; // 1 ghost to create overlap via grown boxes

        // Define a small domain and split into two boxes along x to get neighbor interactions
        IntVect dom_lo(AMREX_D_DECL(0,0,0));
        IntVect dom_hi(AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1));
        Box domain(dom_lo, dom_hi);
        RealBox rb(AMREX_D_DECL(Real(0),Real(0),Real(0)),
                   AMREX_D_DECL(Real(1),Real(1),Real(1)));
        Vector<int> is_per_np(AMREX_SPACEDIM, 0);
        Geometry geom_np(domain, &rb, CoordSys::cartesian, is_per_np.data());
        Vector<int> is_per_p(AMREX_SPACEDIM, 1);
        Geometry geom_p(domain, &rb, CoordSys::cartesian, is_per_p.data());

        // Build an 2x2x2 tiling (8 boxes) to exercise face/edge/corner ghost additions
        IntVect mid = (dom_lo + dom_hi)/2;
        Vector<Box> boxes;
        boxes.reserve(8);
#if (AMREX_SPACEDIM == 3)
        for (int k = 0; k < 2; ++k) {
            for (int j = 0; j < 2; ++j) {
                for (int i = 0; i < 2; ++i) {
                    IntVect lo(AMREX_D_DECL(i==0 ? dom_lo[0] : mid[0]+1,
                                             j==0 ? dom_lo[1] : mid[1]+1,
                                             k==0 ? dom_lo[2] : mid[2]+1));
                    IntVect hi(AMREX_D_DECL(i==0 ? mid[0] : dom_hi[0],
                                             j==0 ? mid[1] : dom_hi[1],
                                             k==0 ? mid[2] : dom_hi[2]));
                    boxes.emplace_back(lo, hi);
                }
            }
        }
#elif (AMREX_SPACEDIM == 2)
        for (int j = 0; j < 2; ++j) {
            for (int i = 0; i < 2; ++i) {
                IntVect lo(AMREX_D_DECL(i==0 ? dom_lo[0] : mid[0]+1,
                                         j==0 ? dom_lo[1] : mid[1]+1,
                                         0));
                IntVect hi(AMREX_D_DECL(i==0 ? mid[0] : dom_hi[0],
                                         j==0 ? mid[1] : dom_hi[1],
                                         0));
                boxes.emplace_back(lo, hi);
            }
        }
#else
        boxes.emplace_back(dom_lo, dom_hi);
#endif
        BoxArray ba(static_cast<size_t>(boxes.size()));
        for (int n = 0; n < static_cast<int>(boxes.size()); ++n) { ba.set(n, boxes[n]); }

        DistributionMapping dm(ba);

        auto run_case = [&](const BoxArray& ba, const DistributionMapping& dm,
                            int ng, const Vector<int>& is_per, const Periodicity& period) {
            // Source values
            MultiFab mf_src(ba, dm, ncomp, ng);
            fill_values(mf_src);
            // Zero out source ghost cells to avoid unintended contributions
            mf_src.setBndry(0.0, 0, ncomp);

            // CPU reference: start from src, zero ghosts, then add neighbor contributions
            MultiFab mf_ref(ba, dm, ncomp, ng);
            amrex::Copy(mf_ref, mf_src, 0, 0, ncomp, ng);
            mf_ref.setBndry(0.0, 0, ncomp);
            IntVect dom_len(AMREX_D_DECL(domain.length(0), domain.length(1), domain.length(2)));

            for (int i = 0; i < ba.size(); ++i) {
                Box v_i = ba.get(i);
                Box gh_i = grow(v_i, ng);
                auto ref_i = mf_ref[i].array();
                for (int j = 0; j < ba.size(); ++j) {
                    Box v_j = ba.get(j);
                    // Iterate periodic image shifts
                    for (int sz = (AMREX_SPACEDIM>=3 && is_per[2]) ? -1 : 0; sz <= ((AMREX_SPACEDIM>=3 && is_per[2]) ? 1 : 0); ++sz) {
                        for (int sy = (AMREX_SPACEDIM>=2 && is_per[1]) ? -1 : 0; sy <= ((AMREX_SPACEDIM>=2 && is_per[1]) ? 1 : 0); ++sy) {
                            for (int sx = is_per[0] ? -1 : 0; sx <= (is_per[0] ? 1 : 0); ++sx) {
                                IntVect shift(AMREX_D_DECL(sx*dom_len[0], sy*dom_len[1], sz*dom_len[2]));
                                if (shift == IntVect::TheZeroVector() && j == i) continue;
                                Box v_js = amrex::shift(v_j, shift);
                                Box region = gh_i & v_js;
                                if (region.ok()) {
                                    auto src_j = mf_src[j].const_array();
                                    int dx = shift[0];
                                    int dy = (AMREX_SPACEDIM>=2 ? shift[1] : 0);
                                    int dz = (AMREX_SPACEDIM>=3 ? shift[2] : 0);
                                    amrex::LoopConcurrentOnCpu(region, ncomp, [=] (int ii, int jj, int kk, int n) noexcept {
                                        ref_i(ii,jj,kk,n) += src_j(ii - dx, jj - dy, kk - dz, n);
                                    });
                                }
                            }
                        }
                    }
                }
            }

            // GPU result using deterministic add
            MultiFab mf_gpu(ba, dm, ncomp, ng);
            amrex::Copy(mf_gpu, mf_src, 0, 0, ncomp, ng);
            mf_gpu.SumBoundaryDeterministic(/*scomp*/0, /*ncomp*/ncomp, IntVect(ng), period);

            // Compare across the full fab box (including ghosts)
            int n_errors = 0;
            for (MFIter mfi(mf_ref); mfi.isValid(); ++mfi) {
                Box const& bx = mf_ref[mfi].box();
                auto const& a = mf_ref[mfi];
                auto const& b = mf_gpu[mfi];
                const int n = 0;
                for (IntVect iv = bx.smallEnd(); iv <= bx.bigEnd(); bx.next(iv)) {
                    if (a(iv,n) != b(iv,n)) {
                        ++n_errors;
                    }
                }
            }
            if (n_errors != 0) {
                amrex::Print() << "ng=" << ng << ": mismatches total=" << n_errors << "\n";
                amrex::Abort("Deterministic add GPU result differs from CPU reference");
            }
        };

        // Non-periodic cases
        run_case(ba, dm, 1, is_per_np, geom_np.periodicity());
        run_case(ba, dm, 2, is_per_np, geom_np.periodicity());
        // Periodic cases (all directions periodic)
        run_case(ba, dm, 1, is_per_p, geom_p.periodicity());
        run_case(ba, dm, 2, is_per_p, geom_p.periodicity());

        // Randomized tilings: generate a few random BoxArray tilings and run the same checks.
        auto make_random_ba = [&](int seed)->BoxArray {
            std::mt19937 gen(seed);
            auto rand_cuts = [&](int len) {
                int kmin = 2, kmax = 4;
                std::uniform_int_distribution<int> dk(kmin, kmax);
                int k = dk(gen);
                std::uniform_int_distribution<int> dpos(0, len-2);
                Vector<int> cuts; cuts.reserve(k-1);
                while ((int)cuts.size() < k-1) {
                    int c = dpos(gen);
                    if (std::find(cuts.begin(), cuts.end(), c) == cuts.end()) cuts.push_back(c);
                }
                std::sort(cuts.begin(), cuts.end());
                Vector<std::pair<int,int>> segs; segs.reserve(k);
                int start = 0;
                for (int idx = 0; idx < (int)cuts.size(); ++idx) {
                    int end = cuts[idx];
                    segs.emplace_back(start, end);
                    start = end+1;
                }
                segs.emplace_back(start, len-1);
                return segs;
            };

            auto segx = rand_cuts(domain.length(0));
#if (AMREX_SPACEDIM >= 2)
            auto segy = rand_cuts(domain.length(1));
#endif
#if (AMREX_SPACEDIM == 3)
            auto segz = rand_cuts(domain.length(2));
#endif
            Vector<Box> rboxes;
            rboxes.reserve(
#if (AMREX_SPACEDIM == 3)
                segx.size()*segy.size()*segz.size()
#elif (AMREX_SPACEDIM == 2)
                segx.size()*segy.size()
#else
                segx.size()
#endif
            );
#if (AMREX_SPACEDIM == 3)
            for (auto zx : segz) {
                for (auto yx : segy) {
                    for (auto xx : segx) {
                        IntVect lo(AMREX_D_DECL(dom_lo[0]+xx.first, dom_lo[1]+yx.first, dom_lo[2]+zx.first));
                        IntVect hi(AMREX_D_DECL(dom_lo[0]+xx.second, dom_lo[1]+yx.second, dom_lo[2]+zx.second));
                        rboxes.emplace_back(lo, hi);
                    }
                }
            }
#elif (AMREX_SPACEDIM == 2)
            for (auto yx : segy) {
                for (auto xx : segx) {
                    IntVect lo(AMREX_D_DECL(dom_lo[0]+xx.first, dom_lo[1]+yx.first, 0));
                    IntVect hi(AMREX_D_DECL(dom_lo[0]+xx.second, dom_lo[1]+yx.second, 0));
                    rboxes.emplace_back(lo, hi);
                }
            }
#else
            for (auto xx : segx) {
                IntVect lo(AMREX_D_DECL(dom_lo[0]+xx.first, 0, 0));
                IntVect hi(AMREX_D_DECL(dom_lo[0]+xx.second, 0, 0));
                rboxes.emplace_back(lo, hi);
            }
#endif
            BoxArray rba(static_cast<size_t>(rboxes.size()));
            for (int n = 0; n < (int)rboxes.size(); ++n) rba.set(n, rboxes[n]);
            return rba;
        };

        for (int trial = 0; trial < 2; ++trial) {
            BoxArray rba = make_random_ba(1234 + trial);
            DistributionMapping rdm(rba);
            // Non-periodic
            run_case(rba, rdm, 1, is_per_np, geom_np.periodicity());
            run_case(rba, rdm, 2, is_per_np, geom_np.periodicity());
            // Periodic (all dirs)
            run_case(rba, rdm, 1, is_per_p, geom_p.periodicity());
            run_case(rba, rdm, 2, is_per_p, geom_p.periodicity());
        }
    }
    amrex::Finalize();
    return 0;
}
