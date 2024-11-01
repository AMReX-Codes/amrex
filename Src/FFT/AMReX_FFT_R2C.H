#ifndef AMREX_FFT_R2C_H_
#define AMREX_FFT_R2C_H_
#include <AMReX_Config.H>

#include <AMReX_MultiFab.H>
#include <AMReX_FFT_Helper.H>
#include <numeric>
#include <tuple>

namespace amrex::FFT
{

/**
 * \brief Discrete Fourier Transform
 *
 * This class supports Fourier transforms between real and complex data. The
 * name R2C indicates that the forward transform converts real data to
 * complex data, while the backward transform converts complex data to real
 * data. It should be noted that both directions of transformation are
 * supported, not just from real to complex. The scaling follows the FFTW
 * convention, where applying the forward transform followed by the backward
 * transform scales the original data by the size of the input array.
 *
 * For more details, we refer the users to
 * https://amrex-codes.github.io/amrex/docs_html/FFT_Chapter.html.
 */
template <typename T = Real, FFT::Direction D = FFT::Direction::both>
class R2C
{
public:
    using MF = std::conditional_t<std::is_same_v<T,Real>,
                                  MultiFab, FabArray<BaseFab<T> > >;
    using cMF = FabArray<BaseFab<GpuComplex<T> > >;

    /**
     * \brief Constructor
     *
     * \param domain the forward domain (i.e., the domain of the real data)
     * \param info optional information
     */
    explicit R2C (Box const& domain, Info const& info = Info{});

    ~R2C ();

    R2C (R2C const&) = delete;
    R2C (R2C &&) = delete;
    R2C& operator= (R2C const&) = delete;
    R2C& operator= (R2C &&) = delete;

    /**
     * \brief Forward and then backward transform
     *
     * This function is available only when this class template is
     * instantiated for transforms in both directions. It's more efficient
     * than calling the forward function that stores the spectral data in a
     * caller provided container followed by the backward function, because
     * this can avoid parallel communication between the internal data and
     * the caller's data container.
     *
     * \param inmf         input data in MultiFab or FabArray<BaseFab<float>>
     * \param outmf        output data in MultiFab or FabArray<BaseFab<float>>
     * \param post_forward a callable object for processing the post-forward
     *                     data before the backward transform. Its interface
     *                     is `(int,int,int,GpuComplex<T>&)`, where the integers
     *                     are indices in the spectral space, and the reference
     *                     to the complex number allows for the modification of
     *                     the spectral data at that location.
     */
    template <typename F, Direction DIR=D,
              std::enable_if_t<DIR == Direction::both, int> = 0>
    void forwardThenBackward (MF const& inmf, MF& outmf, F const& post_forward)
    {
        BL_PROFILE("FFT::R2C::forwardbackward");
        this->forward(inmf);
        this->post_forward_doit(post_forward);
        this->backward(outmf);
    }

    /**
     * \brief Forward transform
     *
     * The output is stored in this object's internal data. This function is
     * not available when this class template is instantiated for
     * backward-only transform.
     *
     * \param inmf input data in MultiFab or FabArray<BaseFab<float>>
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::forward ||
                                                DIR == Direction::both, int> = 0>
    void forward (MF const& inmf);

    /**
     * \brief Forward transform
     *
     * This function is not available when this class template is
     * instantiated for backward-only transform.
     *
     * \param inmf input data in MultiFab or FabArray<BaseFab<float>>
     * \param outmf output data in FabArray<BaseFab<GpuComplex<T>>>
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::forward ||
                                                DIR == Direction::both, int> = 0>
    void forward (MF const& inmf, cMF& outmf);

    /**
     * \brief Backward transform
     *
     * This function is available only when this class template is
     * instantiated for transforms in both directions.
     *
     * \param outmf output data in MultiFab or FabArray<BaseFab<float>>
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::both, int> = 0>
    void backward (MF& outmf);

    /**
     * \brief Backward transform
     *
     * This function is not available when this class template is
     * instantiated for forward-only transform.
     *
     * \param inmf input data in FabArray<BaseFab<GpuComplex<T>>>
     * \param outmf output data in MultiFab or FabArray<BaseFab<float>>
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::backward ||
                                                DIR == Direction::both, int> = 0>
    void backward (cMF const& inmf, MF& outmf);

    /**
     * \brief Get the internal spectral data
     *
     * This function is not available when this class template is
     * instantiated for backward-only transform. For performance reasons,
     * the returned data array does not have the usual ordering of
     * `(x,y,z)`. The order is specified in the second part of the return
     * value.
     */
    template <Direction DIR=D, std::enable_if_t<DIR == Direction::forward ||
                                                DIR == Direction::both, int> = 0>
    std::pair<cMF*,IntVect> getSpectralData ();

    /**
     * \brief Get BoxArray and DistributionMapping for spectral data
     *
     * The returned BoxArray and DistributionMapping can be used to build
     * FabArray<BaseFab<GpuComplex<T>>> for spectral data. The returned
     * BoxArray has the usual order of `(x,y,z)`.
     */
    [[nodiscard]] std::pair<BoxArray,DistributionMapping> getSpectralDataLayout () const;

    // public for cuda
    template <typename F>
    void post_forward_doit (F const& post_forward);

private:

    static std::pair<Plan<T>,Plan<T>> make_c2c_plans (cMF& inout);

    void backward_doit (MF& outmf);

    Plan<T> m_fft_fwd_x{};
    Plan<T> m_fft_bwd_x{};
    Plan<T> m_fft_fwd_y{};
    Plan<T> m_fft_bwd_y{};
    Plan<T> m_fft_fwd_z{};
    Plan<T> m_fft_bwd_z{};

    // Comm meta-data. In the forward phase, we start with (x,y,z),
    // transpose to (y,x,z) and then (z,x,y). In the backward phase, we
    // perform inverse transpose.
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_x2y; // (x,y,z) -> (y,x,z)
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2x; // (y,x,z) -> (x,y,z)
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_y2z; // (y,x,z) -> (z,x,y)
    std::unique_ptr<MultiBlockCommMetaData> m_cmd_z2y; // (z,x,y) -> (y,x,z)
    Swap01 m_dtos_x2y{};
    Swap01 m_dtos_y2x{};
    Swap02 m_dtos_y2z{};
    Swap02 m_dtos_z2y{};

    MF  m_rx;
    cMF m_cx;
    cMF m_cy;
    cMF m_cz;

    std::unique_ptr<char,DataDeleter> m_data_rx_cy;
    std::unique_ptr<char,DataDeleter> m_data_cx_cz;

    Box m_real_domain;
    Box m_spectral_domain_x;
    Box m_spectral_domain_y;
    Box m_spectral_domain_z;

    Info m_info;
};

template <typename T, Direction D>
R2C<T,D>::R2C (Box const& domain, Info const& info)
    : m_real_domain(domain),
      m_spectral_domain_x(IntVect(0), IntVect(AMREX_D_DECL(domain.length(0)/2,
                                                           domain.bigEnd(1),
                                                           domain.bigEnd(2)))),
#if (AMREX_SPACEDIM >= 2)
      m_spectral_domain_y(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(1),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(2)))),
#if (AMREX_SPACEDIM == 3)
      m_spectral_domain_z(IntVect(0), IntVect(AMREX_D_DECL(domain.bigEnd(2),
                                                           domain.length(0)/2,
                                                           domain.bigEnd(1)))),
#endif
#endif
      m_info(info)
{
    BL_PROFILE("FFT::R2C");

    static_assert(std::is_same_v<float,T> || std::is_same_v<double,T>);
    AMREX_ALWAYS_ASSERT(m_real_domain.smallEnd() == 0 &&
                        m_real_domain.length(0) > 1 &&
                        m_real_domain.cellCentered());
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT(m_real_domain.length(2) > 1 || ! m_info.batch_mode);
    AMREX_ALWAYS_ASSERT(m_real_domain.length(1) > 1 || m_real_domain.length(2) == 1);
#else
    AMREX_ALWAYS_ASSERT(! m_info.batch_mode);
#endif

    int myproc = ParallelContext::MyProcSub();
    int nprocs = ParallelContext::NProcsSub();

    //
    // make data containers
    //

    auto bax = amrex::decompose(m_real_domain, nprocs, {AMREX_D_DECL(false,true,true)});
    DistributionMapping dmx = detail::make_iota_distromap(bax.size());
    m_rx.define(bax, dmx, 1, 0, MFInfo().SetAlloc(false));

    {
        BoxList bl = bax.boxList();
        for (auto & b : bl) {
            b.setBig(0, m_spectral_domain_x.bigEnd(0));
        }
        BoxArray cbax(std::move(bl));
        m_cx.define(cbax, dmx, 1, 0, MFInfo().SetAlloc(false));
    }

#if (AMREX_SPACEDIM >= 2)
    DistributionMapping cdmy;
    if (m_real_domain.length(1) > 1) {
        auto cbay = amrex::decompose(m_spectral_domain_y, nprocs, {AMREX_D_DECL(false,true,true)});
        if (cbay.size() == dmx.size()) {
            cdmy = dmx;
        } else {
            cdmy = detail::make_iota_distromap(cbay.size());
        }
        m_cy.define(cbay, cdmy, 1, 0, MFInfo().SetAlloc(false));
    }
#endif

#if (AMREX_SPACEDIM == 3)
    if (m_real_domain.length(1) > 1 &&
        (! m_info.batch_mode && m_real_domain.length(2) > 1))
    {
        auto cbaz = amrex::decompose(m_spectral_domain_z, nprocs, {false,true,true});
        DistributionMapping cdmz;
        if (cbaz.size() == dmx.size()) {
            cdmz = dmx;
        } else if (cbaz.size() == cdmy.size()) {
            cdmz = cdmy;
        } else {
            cdmz = detail::make_iota_distromap(cbaz.size());
        }
         m_cz.define(cbaz, cdmz, 1, 0, MFInfo().SetAlloc(false));
    }
#endif

    m_data_rx_cy = detail::make_mfs_share(m_rx, m_cy);
    m_data_cx_cz = detail::make_mfs_share(m_cx, m_cz);

    //
    // make copiers
    //

#if (AMREX_SPACEDIM >= 2)
    if (! m_cy.empty()) {
        // comm meta-data between x and y phases
        m_cmd_x2y = std::make_unique<MultiBlockCommMetaData>
            (m_cy, m_spectral_domain_y, m_cx, IntVect(0), m_dtos_x2y);
        m_cmd_y2x = std::make_unique<MultiBlockCommMetaData>
            (m_cx, m_spectral_domain_x, m_cy, IntVect(0), m_dtos_y2x);
    }
#endif
#if (AMREX_SPACEDIM == 3)
    if (! m_cz.empty()) {
        // comm meta-data between y and z phases
        m_cmd_y2z = std::make_unique<MultiBlockCommMetaData>
            (m_cz, m_spectral_domain_z, m_cy, IntVect(0), m_dtos_y2z);
        m_cmd_z2y = std::make_unique<MultiBlockCommMetaData>
            (m_cy, m_spectral_domain_y, m_cz, IntVect(0), m_dtos_z2y);
    }
#endif

    //
    // make plans
    //

    if (myproc < m_rx.size())
    {
        Box const& box = m_rx.box(myproc);
        auto* pr = m_rx[myproc].dataPtr();
        auto* pc = (typename Plan<T>::VendorComplex *)m_cx[myproc].dataPtr();
#ifdef AMREX_USE_SYCL
        m_fft_fwd_x.template init_r2c<Direction::forward>(box, pr, pc);
        m_fft_bwd_x = m_fft_fwd_x;
#else
        if constexpr (D == Direction::both || D == Direction::forward) {
            m_fft_fwd_x.template init_r2c<Direction::forward>(box, pr, pc);
        }
        if constexpr (D == Direction::both || D == Direction::backward) {
            m_fft_bwd_x.template init_r2c<Direction::backward>(box, pr, pc);
        }
#endif
    }

#if (AMREX_SPACEDIM >= 2)
    if (! m_cy.empty()) {
        std::tie(m_fft_fwd_y, m_fft_bwd_y) = make_c2c_plans(m_cy);
    }
#endif
#if (AMREX_SPACEDIM == 3)
    if (! m_cz.empty()) {
        std::tie(m_fft_fwd_z, m_fft_bwd_z) = make_c2c_plans(m_cz);
    }
#endif
}

template <typename T, Direction D>
R2C<T,D>::~R2C<T,D> ()
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

template <typename T, Direction D>
template <Direction DIR, std::enable_if_t<DIR == Direction::forward ||
                                          DIR == Direction::both, int> >
void R2C<T,D>::forward (MF const& inmf)
{
    BL_PROFILE("FFT::R2C::forward(in)");

    m_rx.ParallelCopy(inmf, 0, 0, 1);
    m_fft_fwd_x.template compute_r2c<Direction::forward>();

    if (                          m_cmd_x2y) {
        ParallelCopy(m_cy, m_cx, *m_cmd_x2y, 0, 0, 1, m_dtos_x2y);
    }
    m_fft_fwd_y.template compute_c2c<Direction::forward>();

    if (                          m_cmd_y2z) {
        ParallelCopy(m_cz, m_cy, *m_cmd_y2z, 0, 0, 1, m_dtos_y2z);
    }
    m_fft_fwd_z.template compute_c2c<Direction::forward>();
}

template <typename T, Direction D>
template <Direction DIR, std::enable_if_t<DIR == Direction::both, int> >
void R2C<T,D>::backward (MF& outmf)
{
    backward_doit(outmf);
}

template <typename T, Direction D>
void R2C<T,D>::backward_doit (MF& outmf)
{
    BL_PROFILE("FFT::R2C::backward(out)");

    m_fft_bwd_z.template compute_c2c<Direction::backward>();
    if (                          m_cmd_z2y) {
        ParallelCopy(m_cy, m_cz, *m_cmd_z2y, 0, 0, 1, m_dtos_z2y);
    }

    m_fft_bwd_y.template compute_c2c<Direction::backward>();
    if (                          m_cmd_y2x) {
        ParallelCopy(m_cx, m_cy, *m_cmd_y2x, 0, 0, 1, m_dtos_y2x);
    }

    m_fft_bwd_x.template compute_r2c<Direction::backward>();
    outmf.ParallelCopy(m_rx, 0, 0, 1);
}

template <typename T, Direction D>
std::pair<Plan<T>, Plan<T>>
R2C<T,D>::make_c2c_plans (cMF& inout)
{
    Plan<T> fwd;
    Plan<T> bwd;

    auto* fab = detail::get_fab(inout);
    if (!fab) { return {fwd, bwd};}

    Box const& box = fab->box();
    auto* pio = (typename Plan<T>::VendorComplex *)fab->dataPtr();

#ifdef AMREX_USE_SYCL
    fwd.template init_c2c<Direction::forward>(box, pio);
    bwd = fwd;
#else
    if constexpr (D == Direction::both || D == Direction::forward) {
        fwd.template init_c2c<Direction::forward>(box, pio);
    }
    if constexpr (D == Direction::both || D == Direction::backward) {
        bwd.template init_c2c<Direction::backward>(box, pio);
    }
#endif

    return {fwd, bwd};
}

template <typename T, Direction D>
template <typename F>
void R2C<T,D>::post_forward_doit (F const& post_forward)
{
    if (m_info.batch_mode) {
        amrex::Abort("xxxxx todo: post_forward");
    } else {
        if (                           ! m_cz.empty()) {
            auto* spectral_fab = detail::get_fab(m_cz);
            if (spectral_fab) {
                auto const& a = spectral_fab->array(); // m_cz's ordering is z,x,y
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int iz, int jx, int ky)
                {
                    post_forward(jx,ky,iz,a(iz,jx,ky));
                });
            }
        } else if (                    ! m_cy.empty()) {
            auto* spectral_fab = detail::get_fab(m_cy);
            if (spectral_fab) {
                auto const& a = spectral_fab->array(); // m_cy's ordering is y,x,z
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int iy, int jx, int k)
                {
                    post_forward(jx,iy,k,a(iy,jx,k));
                });
            }
        } else {
            auto* spectral_fab = detail::get_fab(m_cx);
            if (spectral_fab) {
                auto const& a = spectral_fab->array();
                ParallelFor(spectral_fab->box(),
                [=] AMREX_GPU_DEVICE (int i, int j, int k)
                {
                    post_forward(i,j,k,a(i,j,k));
                });
            }
        }
    }
}

template <typename T, Direction D>
template <Direction DIR, std::enable_if_t<DIR == Direction::forward ||
                                          DIR == Direction::both, int> >
std::pair<typename R2C<T,D>::cMF *, IntVect>
R2C<T,D>::getSpectralData ()
{
    if (!m_cz.empty()) {
        return std::make_pair(&m_cz, IntVect{AMREX_D_DECL(2,0,1)});
    } else if (!m_cy.empty()) {
        return std::make_pair(&m_cy, IntVect{AMREX_D_DECL(1,0,2)});
    } else {
        return std::make_pair(&m_cx, IntVect{AMREX_D_DECL(0,1,2)});
    }
}

template <typename T, Direction D>
template <Direction DIR, std::enable_if_t<DIR == Direction::forward ||
                                          DIR == Direction::both, int> >
void R2C<T,D>::forward (MF const& inmf, cMF& outmf)
{
    BL_PROFILE("FFT::R2C::forward(inout)");

    forward(inmf);
    if (!m_cz.empty()) { // m_cz's order (z,x,y) -> (x,y,z)
        RotateBwd dtos{};
        MultiBlockCommMetaData cmd
            (outmf, m_spectral_domain_x, m_cz, IntVect(0), dtos);
        ParallelCopy(outmf, m_cz, cmd, 0, 0, 1, dtos);
    } else if (!m_cy.empty()) { // m_cy's order (y,x,z) -> (x,y,z)
        MultiBlockCommMetaData cmd
            (outmf, m_spectral_domain_x, m_cy, IntVect(0), m_dtos_y2x);
        ParallelCopy(outmf, m_cy, cmd, 0, 0, 1, m_dtos_y2x);
    } else {
        outmf.ParallelCopy(m_cx, 0, 0, 1);
    }
}

template <typename T, Direction D>
template <Direction DIR, std::enable_if_t<DIR == Direction::backward ||
                                          DIR == Direction::both, int> >
void R2C<T,D>::backward (cMF const& inmf, MF& outmf)
{
    BL_PROFILE("FFT::R2C::backward(inout)");

    if (!m_cz.empty()) { // (x,y,z) -> m_cz's order (z,x,y)
        RotateFwd dtos{};
        MultiBlockCommMetaData cmd
            (m_cz, m_spectral_domain_z, inmf, IntVect(0), dtos);
        ParallelCopy(m_cz, inmf, cmd, 0, 0, 1, dtos);
    } else if (!m_cy.empty()) { // (x,y,z) -> m_cy's ordering (y,x,z)
        MultiBlockCommMetaData cmd
            (m_cy, m_spectral_domain_y, inmf, IntVect(0), m_dtos_x2y);
        ParallelCopy(m_cy, inmf, cmd, 0, 0, 1, m_dtos_x2y);
    } else {
        m_cx.ParallelCopy(inmf, 0, 0, 1);
    }
    backward_doit(outmf);
}

template <typename T, Direction D>
std::pair<BoxArray,DistributionMapping>
R2C<T,D>::getSpectralDataLayout () const
{
#if (AMREX_SPACEDIM == 3)
    if (!m_cz.empty()) {
        BoxList bl = m_cz.boxArray().boxList();
        for (auto& b : bl) {
            auto lo = b.smallEnd();
            auto hi = b.bigEnd();
            std::swap(lo[0], lo[1]);
            std::swap(lo[1], lo[2]);
            std::swap(hi[0], hi[1]);
            std::swap(hi[1], hi[2]);
            b.setSmall(lo);
            b.setBig(hi);
        }
        return std::make_pair(BoxArray(std::move(bl)), m_cz.DistributionMap());
    } else
#endif
#if (AMREX_SPACEDIM >= 2)
    if (!!m_cy.empty()) {
        BoxList bl = m_cy.boxArray().boxList();
        for (auto& b : bl) {
            auto lo = b.smallEnd();
            auto hi = b.bigEnd();
            std::swap(lo[0], lo[1]);
            std::swap(hi[0], hi[1]);
            b.setSmall(lo);
            b.setBig(hi);
        }
        return std::make_pair(BoxArray(std::move(bl)), m_cy.DistributionMap());
    } else
#endif
    {
        return std::make_pair(m_cx.boxArray(), m_cx.DistributionMap());
    }
}

}

#endif
