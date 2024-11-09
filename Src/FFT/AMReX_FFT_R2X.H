#ifndef AMREX_FFT_R2X_H_
#define AMREX_FFT_R2X_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_FFT_Helper.H>
#include <algorithm>
#include <numeric>
#include <tuple>

namespace amrex::FFT
{

/**
 * \brief Discrete Fourier Transform
 *
 * This class supports Fourier transforms including cosine and sine
 * transforms.
 */
template <typename T = Real>
class R2X
{
public:
    using MF = std::conditional_t<std::is_same_v<T,Real>,
                                  MultiFab, FabArray<BaseFab<T> > >;
    using cMF = FabArray<BaseFab<GpuComplex<T> > >;

    R2X (Box const& domain,
         Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc,
         Info const& info = Info{});

    ~R2X ();

    R2X (R2X const&) = delete;
    R2X (R2X &&) = delete;
    R2X& operator= (R2X const&) = delete;
    R2X& operator= (R2X &&) = delete;

    [[nodiscard]] T scalingFactor () const;

    template <typename F>
    void forwardThenBackward (MF const& inmf, MF& outmf, F const& post_forward);

    // public for cuda
    template <int dim, typename FAB, typename F>
    void post_forward_doit (FAB* fab, F const& f);

private:
    Box m_dom_0;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;

    Plan<T> m_fft_fwd_x{};
    Plan<T> m_fft_bwd_x{};
    Plan<T> m_fft_fwd_y{};
    Plan<T> m_fft_bwd_y{};
    Plan<T> m_fft_fwd_z{};
    Plan<T> m_fft_bwd_z{};

    std::unique_ptr<MultiBlockCommMetaData> m_cmd_cx2cy;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_rx2ry;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_cy2cz;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_ry2rz;

    std::unique_ptr<MultiBlockCommMetaData> m_cmd_cy2cx;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_ry2rx;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_cz2cy;
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_rz2ry;

    Swap01 m_dtos_x2y{};
    Swap01 m_dtos_y2x{};
    Swap02 m_dtos_y2z{};
    Swap02 m_dtos_z2y{};

    MF m_rx;
    MF m_ry;
    MF m_rz;
    cMF m_cx;
    cMF m_cy;
    cMF m_cz;

    std::unique_ptr<char,DataDeleter> m_data_1;
    std::unique_ptr<char,DataDeleter> m_data_2;

    Box m_dom_rx;
    Box m_dom_ry;
    Box m_dom_rz;
    Box m_dom_cx;
    Box m_dom_cy;
    Box m_dom_cz;

    Info m_info;
};

template <typename T>
R2X<T>::R2X (Box const& domain,
             Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc,
             Info const& info)
    : m_dom_0(domain),
      m_bc(bc),
      m_info(info)
{
    BL_PROFILE("FFT::R2X");

    static_assert(std::is_same_v<float,T> || std::is_same_v<double,T>);
    AMREX_ALWAYS_ASSERT(domain.smallEnd() == 0 &&
                        domain.length(0) > 1 &&
                        domain.cellCentered());
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT(domain.length(1) > 1 || domain.length(2) == 1);
#endif
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT(domain.length(idim) > 1);
        if (bc[idim].first == Boundary::periodic ||
            bc[idim].second == Boundary::periodic) {
            AMREX_ALWAYS_ASSERT(bc[idim].first == bc[idim].second);
        }
    }

    int myproc = ParallelContext::MyProcSub();
    int nprocs = std::min(ParallelContext::NProcsSub(), m_info.nprocs);

    //
    // make data containers
    //

    m_dom_rx = m_dom_0;
    auto bax = amrex::decompose(m_dom_rx, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping dmx = detail::make_iota_distromap(bax.size());
    m_rx.define(bax, dmx, 1, 0, MFInfo().SetAlloc(false));

    // x-direction
    if (bc[0].first == Boundary::periodic) {
        // x-fft: r2c(m_rx->m_cx)
        m_dom_cx = Box(IntVect(0), IntVect(AMREX_D_DECL(domain.length(0)/2,
                                                        domain.bigEnd(1),
                                                        domain.bigEnd(2))));
        BoxList bl = bax.boxList();
        for (auto & b : bl) {
            b.setBig(0, m_dom_cx.bigEnd(0));
        }
        BoxArray cbax(std::move(bl));
        m_cx.define(cbax, dmx, 1, 0, MFInfo().SetAlloc(false));
    } // else: x-fft: r2r(m_rx)

#if (AMREX_SPACEDIM >= 2)
    if (domain.length(1) > 1) {
        if (! m_cx.empty()) {
            // copy(m_cx->m_cy)
            m_dom_cy = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_cx.bigEnd(1),
                                                            m_dom_cx.bigEnd(0),
                                                            m_dom_cx.bigEnd(2))));
            auto ba = amrex::decompose(m_dom_cy, nprocs, {AMREX_D_DECL(false,true,true)});
            DistributionMapping dm;
            if (ba.size() == m_cx.size()) {
                dm = m_cx.DistributionMap();
            } else {
                dm = detail::make_iota_distromap(ba.size());
            }
            m_cy.define(ba, dm, 1, 0, MFInfo().SetAlloc(false));
            // if bc[1] is periodic:
            //     c2c(m_cy->m_cy)
            // else:
            //     r2r(m_cy.re) & r2r(m_cy.im)
        } else {
            // copy(m_rx->m_ry)
            m_dom_ry = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_rx.bigEnd(1),
                                                            m_dom_rx.bigEnd(0),
                                                            m_dom_rx.bigEnd(2))));
            auto ba = amrex::decompose(m_dom_ry, nprocs, {AMREX_D_DECL(false,true,true)});
            DistributionMapping dm;
            if (ba.size() == m_rx.size()) {
                dm = m_rx.DistributionMap();
            } else {
                dm = detail::make_iota_distromap(ba.size());
            }
            m_ry.define(ba, dm, 1, 0, MFInfo().SetAlloc(false));
            // if bc[1] is periodic:
            //     r2c(m_ry->m_cy)
            // else:
            //     r2r(m_ry)
            if (bc[1].first == Boundary::periodic) {
                m_dom_cy = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_ry.length(0)/2,
                                                                m_dom_ry.bigEnd(1),
                                                                m_dom_ry.bigEnd(2))));
                BoxList bl = ba.boxList();
                for (auto & b : bl) {
                    b.setBig(0, m_dom_cy.bigEnd(0));
                }
                BoxArray cba(std::move(bl));
                m_cy.define(cba, dm, 1, 0, MFInfo().SetAlloc(false));
            }
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    if (domain.length(2) > 1) {
        if (! m_cy.empty()) {
            // copy(m_cy, m_cz)
            m_dom_cz = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_cy.bigEnd(2),
                                                            m_dom_cy.bigEnd(1),
                                                            m_dom_cy.bigEnd(0))));
            auto ba = amrex::decompose(m_dom_cz, nprocs, {AMREX_D_DECL(false,true,true)});
            DistributionMapping dm;
            if (ba.size() == m_cy.size()) {
                dm = m_cy.DistributionMap();
            } else {
                dm = detail::make_iota_distromap(ba.size());
            }
            m_cz.define(ba, dm, 1, 0, MFInfo().SetAlloc(false));
            // if bc[2] is periodic:
            //     c2c(m_cz->m_cz)
            // else:
            //     r2r(m_cz.re) & r2r(m_cz.im)
        } else {
            // copy(m_ry, m_rz)
            m_dom_rz = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_ry.bigEnd(2),
                                                            m_dom_ry.bigEnd(1),
                                                            m_dom_ry.bigEnd(0))));
            auto ba = amrex::decompose(m_dom_rz, nprocs, {AMREX_D_DECL(false,true,true)});
            DistributionMapping dm;
            if (ba.size() == m_ry.size()) {
                dm = m_ry.DistributionMap();
            } else {
                dm = detail::make_iota_distromap(ba.size());
            }
            m_rz.define(ba, dm, 1, 0, MFInfo().SetAlloc(false));
            // if bc[2] is periodic:
            //     r2c(m_rz->m_cz)
            // else:
            //     r2r(m_rz)
            if (bc[2].first == Boundary::periodic) {
                m_dom_cz = Box(IntVect(0), IntVect(AMREX_D_DECL(m_dom_rz.length(0)/2,
                                                                m_dom_rz.bigEnd(1),
                                                                m_dom_rz.bigEnd(2))));
                BoxList bl = ba.boxList();
                for (auto & b : bl) {
                    b.setBig(0, m_dom_cz.bigEnd(0));
                }
                BoxArray cba(std::move(bl));
                m_cz.define(cba, dm, 1, 0, MFInfo().SetAlloc(false));
            }
        }
    }
#endif

    // There are several different execution paths.
    //
    // (1) x-r2c(m_rx->m_cx), copy(m_cx->m_cy), y-fft(m_cy),
    //     copy(m_cy->m_cz), z-fft(m_cz)
    //     In this case, we have m_rx, m_cx, m_cy, & m_cz.
    //     we can alias(m_rx,m_cy) and alias(m_cx,m_cz).
    //
    // (2) x_r2r(m_rx), copy(m_rx->m_ry), y-r2c(m_ry->m_cy),
    //     copy(m_cy->m_cz), z-fft(m_cz)
    //     In this case, we have m_rx, m_ry, m_cy, & m_cz.
    //     We can alias(m_rx,m_cy) and alias(m_ry,m_cz).
    //
    // (3) x_r2r(m_rx), copy(m_rx->m_ry), y-r2r(m_ry),
    //     copy(m_ry->m_rz), z-r2c(m_rz->m_rz)
    //     In this case, we have m_rx, m_ry, m_rz, & m_cz
    //     We can alias(m_rx,m_rz) and alias(m_ry,m_cz)
    //
    // (4) x_r2r(m_rx), copy(m_rx->m_ry), y-r2r(m_ry),
    //     copy(m_ry->m_rz), z-r2r(m_rz)
    //     In this case, we have m_rx, m_ry, & m_rz.
    //     We can alias(m_rx,m_rz).

    if (! m_cx.empty()) {
        m_data_1 = detail::make_mfs_share(m_rx, m_cy);
        m_data_2 = detail::make_mfs_share(m_cx, m_cz);
    } else if (! m_cy.empty()) {
        m_data_1 = detail::make_mfs_share(m_rx, m_cy);
        m_data_2 = detail::make_mfs_share(m_ry, m_cz);
    } else if (! m_cz.empty()) {
        m_data_1 = detail::make_mfs_share(m_rx, m_rz);
        m_data_2 = detail::make_mfs_share(m_ry, m_cz);
    } else {
        m_data_1 = detail::make_mfs_share(m_rx, m_rz);
        m_data_2 = detail::make_mfs_share(m_ry, m_cz); // It's okay m_cz is empty.
    }

    //
    // make copiers
    //

#if (AMREX_SPACEDIM >= 2)
    if (domain.length(1) > 1) {
        if (! m_cx.empty()) {
            // copy(m_cx->m_cy)
            m_cmd_cx2cy = std::make_unique<MultiBlockCommMetaData>
                (m_cy, m_dom_cy, m_cx, IntVect(0), m_dtos_x2y);
            m_cmd_cy2cx = std::make_unique<MultiBlockCommMetaData>
                (m_cx, m_dom_cx, m_cy, IntVect(0), m_dtos_y2x);
        } else {
            // copy(m_rx->m_ry)
            m_cmd_rx2ry = std::make_unique<MultiBlockCommMetaData>
                (m_ry, m_dom_ry, m_rx, IntVect(0), m_dtos_x2y);
            m_cmd_ry2rx = std::make_unique<MultiBlockCommMetaData>
                (m_rx, m_dom_rx, m_ry, IntVect(0), m_dtos_y2x);
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    if (domain.length(2) > 1) {
        if (! m_cy.empty()) {
            // copy(m_cy, m_cz)
            m_cmd_cy2cz = std::make_unique<MultiBlockCommMetaData>
                (m_cz, m_dom_cz, m_cy, IntVect(0), m_dtos_y2z);
            m_cmd_cz2cy = std::make_unique<MultiBlockCommMetaData>
                (m_cy, m_dom_cy, m_cz, IntVect(0), m_dtos_z2y);
        } else {
            // copy(m_ry, m_rz)
            m_cmd_ry2rz = std::make_unique<MultiBlockCommMetaData>
                (m_rz, m_dom_rz, m_ry, IntVect(0), m_dtos_y2z);
            m_cmd_rz2ry = std::make_unique<MultiBlockCommMetaData>
                (m_ry, m_dom_ry, m_rz, IntVect(0), m_dtos_z2y);
        }
    }
#endif

    //
    // make plans
    //

    using VendorComplex = typename Plan<T>::VendorComplex;

    if (myproc < m_rx.size())
    {
        Box const& box = m_rx.box(myproc);
        auto* pf = m_rx[myproc].dataPtr();
        if (bc[0].first == Boundary::periodic) {
            auto* pb = (VendorComplex*) m_cx[myproc].dataPtr();
            m_fft_fwd_x.template init_r2c<Direction::forward>(box, pf, pb);
#if defined(AMREX_USE_SYCL)
            m_fft_bwd_x = m_fft_fwd_x;
#else
            m_fft_bwd_x.template init_r2c<Direction::backward>(box, pf, pb);
#endif
        } else {
            m_fft_fwd_x.template init_r2r<Direction::forward>(box, pf, bc[0]);
#if defined(AMREX_USE_GPU)
            if ((bc[0].first == Boundary::even && bc[0].second == Boundary::odd) ||
                (bc[0].first == Boundary::odd && bc[0].second == Boundary::even)) {
                m_fft_bwd_x = m_fft_fwd_x;
            } else
#endif
            {
                m_fft_bwd_x.template init_r2r<Direction::backward>(box, pf, bc[0]);
            }
        }
    }

#if (AMREX_SPACEDIM >= 2)
    if (m_ry.empty() && m_bc[1].first == Boundary::periodic) {
        if (myproc < m_cy.size()) {
            Box const& box = m_cy.box(myproc);
            auto* p = (VendorComplex *)m_cy[myproc].dataPtr();
            m_fft_fwd_y.template init_c2c<Direction::forward>(box, p);
#if defined(AMREX_USE_SYCL)
            m_fft_bwd_y = m_fft_fwd_y;
#else
            m_fft_bwd_y.template init_c2c<Direction::backward>(box, p);
#endif
        }
    } else if (!m_ry.empty() && m_bc[1].first == Boundary::periodic) {
        if (myproc < m_ry.size()) {
            Box const& box = m_ry.box(myproc);
            auto* pr =                 m_ry[myproc].dataPtr();
            auto* pc = (VendorComplex*)m_cy[myproc].dataPtr();
            m_fft_fwd_y.template init_r2c<Direction::forward>(box, pr, pc);
#if defined(AMREX_USE_SYCL)
            m_fft_bwd_y = m_fft_fwd_y;
#else
            m_fft_bwd_y.template init_r2c<Direction::backward>(box, pr, pc);
#endif
        }
    } else if (!m_cy.empty()) {
        if (myproc < m_cy.size()) {
            Box const& box = m_cy.box(myproc);
            auto* p = (VendorComplex*) m_cy[myproc].dataPtr();
            m_fft_fwd_y.template init_r2r<Direction::forward>(box, p, bc[1]);
#if defined(AMREX_USE_GPU)
            if ((bc[1].first == Boundary::even && bc[1].second == Boundary::odd) ||
                (bc[1].first == Boundary::odd && bc[1].second == Boundary::even)) {
                m_fft_bwd_y = m_fft_fwd_y;
            } else
#endif
            {
                m_fft_bwd_y.template init_r2r<Direction::backward>(box, p, bc[1]);
            }
        }
    } else {
        if (myproc < m_ry.size()) {
            Box const& box = m_ry.box(myproc);
            auto* p = m_ry[myproc].dataPtr();
            m_fft_fwd_y.template init_r2r<Direction::forward>(box, p, bc[1]);
#if defined(AMREX_USE_GPU)
            if ((bc[1].first == Boundary::even && bc[1].second == Boundary::odd) ||
                (bc[1].first == Boundary::odd && bc[1].second == Boundary::even)) {
                m_fft_bwd_y = m_fft_fwd_y;
            } else
#endif
            {
                m_fft_bwd_y.template init_r2r<Direction::backward>(box, p, bc[1]);
            }
        }
    }
#endif

#if (AMREX_SPACEDIM == 3)
    if (m_rz.empty() && m_bc[2].first == Boundary::periodic) {
        if (myproc < m_cz.size()) {
            Box const& box = m_cz.box(myproc);
            auto* p = (VendorComplex*)m_cz[myproc].dataPtr();
            m_fft_fwd_z.template init_c2c<Direction::forward>(box, p);
#if defined(AMREX_USE_SYCL)
            m_fft_bwd_z = m_fft_fwd_z;
#else
            m_fft_bwd_z.template init_c2c<Direction::backward>(box, p);
#endif
        }
    } else if (!m_rz.empty() && m_bc[2].first == Boundary::periodic) {
        if (myproc < m_rz.size()) {
            Box const& box = m_rz.box(myproc);
            auto* pr =                 m_rz[myproc].dataPtr();
            auto* pc = (VendorComplex*)m_cz[myproc].dataPtr();
            m_fft_fwd_z.template init_r2c<Direction::forward>(box, pr, pc);
#if defined(AMREX_USE_SYCL)
            m_fft_bwd_z = m_fft_fwd_z;
#else
            m_fft_bwd_z.template init_r2c<Direction::backward>(box, pr, pc);
#endif
        }
    } else if (!m_cz.empty()) {
        if (myproc < m_cz.size()) {
            Box const& box = m_cz.box(myproc);
            auto* p = (VendorComplex*) m_cz[myproc].dataPtr();
            m_fft_fwd_z.template init_r2r<Direction::forward>(box, p, bc[2]);
#if defined(AMREX_USE_GPU)
            if ((bc[2].first == Boundary::even && bc[2].second == Boundary::odd) ||
                (bc[2].first == Boundary::odd && bc[2].second == Boundary::even)) {
                m_fft_bwd_z = m_fft_fwd_z;
            } else
#endif
            {
                m_fft_bwd_z.template init_r2r<Direction::backward>(box, p, bc[2]);
            }
        }
    } else {
        if (myproc < m_rz.size()) {
            Box const& box = m_rz.box(myproc);
            auto* p = m_rz[myproc].dataPtr();
            m_fft_fwd_z.template init_r2r<Direction::forward>(box, p, bc[2]);
#if defined(AMREX_USE_GPU)
            if ((bc[2].first == Boundary::even && bc[2].second == Boundary::odd) ||
                (bc[2].first == Boundary::odd && bc[2].second == Boundary::even)) {
                m_fft_bwd_z = m_fft_fwd_z;
            } else
#endif
            {
                m_fft_bwd_z.template init_r2r<Direction::backward>(box, p, bc[2]);
            }
        }
    }
#endif
}

template <typename T>
R2X<T>::~R2X ()
{
    if (m_fft_bwd_x.plan != m_fft_fwd_x.plan) {
        m_fft_bwd_x.destroy();
    }
    if (m_fft_bwd_y.plan != m_fft_fwd_y.plan) {
        m_fft_bwd_y.destroy();
    }
    if (m_fft_bwd_z.plan != m_fft_fwd_z.plan) {
        m_fft_bwd_z.destroy();
    }
    m_fft_fwd_x.destroy();
    m_fft_fwd_y.destroy();
    m_fft_fwd_z.destroy();
}

template <typename T>
T R2X<T>::scalingFactor () const
{
    auto r = m_dom_0.numPts();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        if (m_bc[idim].first != Boundary::periodic
            && m_dom_0.length(idim) > 1)
        {
            r *= 2;
        }
    }
    return T(1)/T(r);
}

template <typename T>
template <typename F>
void R2X<T>::forwardThenBackward (MF const& inmf, MF& outmf, F const& post_forward)
{
    BL_PROFILE("FFT::R2X::forwardbackward");

    // forward

    m_rx.ParallelCopy(inmf, 0, 0, 1);
    if (m_bc[0].first == Boundary::periodic) {
        m_fft_fwd_x.template compute_r2c<Direction::forward>();
    } else {
        m_fft_fwd_x.template compute_r2r<Direction::forward>();
    }

#if (AMREX_SPACEDIM >= 2)
    if (                          m_cmd_cx2cy) {
        ParallelCopy(m_cy, m_cx, *m_cmd_cx2cy, 0, 0, 1, m_dtos_x2y);
    } else if (                   m_cmd_rx2ry) {
        ParallelCopy(m_ry, m_rx, *m_cmd_rx2ry, 0, 0, 1, m_dtos_x2y);
    }
    if (m_bc[1].first != Boundary::periodic)
    {
        m_fft_fwd_y.template compute_r2r<Direction::forward>();
    }
    else if (m_bc[0].first == Boundary::periodic)
    {
        m_fft_fwd_y.template compute_c2c<Direction::forward>();
    }
    else
    {
        m_fft_fwd_y.template compute_r2c<Direction::forward>();
    }
#endif

#if (AMREX_SPACEDIM == 3)
    if (                          m_cmd_cy2cz) {
        ParallelCopy(m_cz, m_cy, *m_cmd_cy2cz, 0, 0, 1, m_dtos_y2z);
    } else if (                   m_cmd_ry2rz) {
        ParallelCopy(m_rz, m_ry, *m_cmd_ry2rz, 0, 0, 1, m_dtos_y2z);
    }
    if (m_bc[2].first != Boundary::periodic)
    {
        m_fft_fwd_z.template compute_r2r<Direction::forward>();
    }
    else if (m_bc[0].first == Boundary::periodic ||
             m_bc[1].first == Boundary::periodic)
    {
        m_fft_fwd_z.template compute_c2c<Direction::forward>();
    }
    else
    {
        m_fft_fwd_z.template compute_r2c<Direction::forward>();
    }
#endif

    // post-forward

    int actual_dim = AMREX_SPACEDIM;
#if (AMREX_SPACEDIM >= 2)
    if (m_dom_0.length(1) == 1) { actual_dim = 1; }
#endif
#if (AMREX_SPACEDIM == 3)
    if ((m_dom_0.length(2) == 1) && (m_dom_0.length(1) > 1)) { actual_dim = 2; }
#endif

    if (actual_dim == 1) {
        if (m_cx.empty()) {
            post_forward_doit<0>(detail::get_fab(m_rx), post_forward);
        } else {
            post_forward_doit<0>(detail::get_fab(m_cx), post_forward);
        }
    }
#if (AMREX_SPACEDIM >= 2)
    else if (actual_dim == 2) {
        if (m_cy.empty()) {
            post_forward_doit<1>(detail::get_fab(m_ry), post_forward);
        } else {
            post_forward_doit<1>(detail::get_fab(m_cy), post_forward);
        }
    }
#endif
#if (AMREX_SPACEDIM == 3)
    else if (actual_dim == 3) {
        if (m_cz.empty()) {
            post_forward_doit<2>(detail::get_fab(m_rz), post_forward);
        } else {
            post_forward_doit<2>(detail::get_fab(m_cz), post_forward);
        }
    }
#endif

    // backward

#if (AMREX_SPACEDIM == 3)
    if (m_bc[2].first != Boundary::periodic)
    {
        m_fft_bwd_z.template compute_r2r<Direction::backward>();
    }
    else if (m_bc[0].first == Boundary::periodic ||
             m_bc[1].first == Boundary::periodic)
    {
        m_fft_bwd_z.template compute_c2c<Direction::backward>();
    }
    else
    {
        m_fft_bwd_z.template compute_r2c<Direction::backward>();
    }
    if (                          m_cmd_cz2cy) {
        ParallelCopy(m_cy, m_cz, *m_cmd_cz2cy, 0, 0, 1, m_dtos_z2y);
    } else if (                   m_cmd_rz2ry) {
        ParallelCopy(m_ry, m_rz, *m_cmd_rz2ry, 0, 0, 1, m_dtos_z2y);
    }
#endif

#if (AMREX_SPACEDIM >= 2)
    if (m_bc[1].first != Boundary::periodic)
    {
        m_fft_bwd_y.template compute_r2r<Direction::backward>();
    }
    else if (m_bc[0].first == Boundary::periodic)
    {
        m_fft_bwd_y.template compute_c2c<Direction::backward>();
    }
    else
    {
        m_fft_bwd_y.template compute_r2c<Direction::backward>();
    }
    if (                          m_cmd_cy2cx) {
        ParallelCopy(m_cx, m_cy, *m_cmd_cy2cx, 0, 0, 1, m_dtos_y2x);
    } else if (                   m_cmd_ry2rx) {
        ParallelCopy(m_rx, m_ry, *m_cmd_ry2rx, 0, 0, 1, m_dtos_y2x);
    }
#endif

    if (m_bc[0].first == Boundary::periodic) {
        m_fft_bwd_x.template compute_r2c<Direction::backward>();
    } else {
        m_fft_bwd_x.template compute_r2r<Direction::backward>();
    }
    outmf.ParallelCopy(m_rx, 0, 0, 1);
}

template <typename T>
template <int dim, typename FAB, typename F>
void R2X<T>::post_forward_doit (FAB* fab, F const& f)
{
    if (fab) {
        auto const& a = fab->array();
        ParallelFor(fab->box(),
        [f=f,a=a] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if constexpr (dim == 0) {
                f(i,j,k,a(i,j,k));
            } else if constexpr (dim == 1) {
                f(j,i,k,a(i,j,k));
            } else {
                f(j,k,i,a(i,j,k));
            }
        });
    }
}

}

#endif
