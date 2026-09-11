#ifndef AMREX_FFT_OPENBC_SOLVER_H_
#define AMREX_FFT_OPENBC_SOLVER_H_

#include <AMReX_FFT_R2C.H>
#include <AMReX_FFT_R2X.H>

namespace amrex::FFT
{

/**
 * \file AMReX_FFT_OpenBCSolver.H
 *
 * FFT-based solver for open-boundary convolution problems via Green's functions.
 */

/**
 * \ingroup amrex_fft
 * \brief Convolution-based solver for open boundary conditions using Green's functions.
 *
 * The class wraps an internal `R2C` plan and multiplies spectral data by a user-provided
 * Green's function to solve Laplace/Poisson problems on a rectangular domain.
 * It currently solves one right-hand-side component per call; the supplied
 * `FFT::Info` must have `batch_size == 1`.
 */
template <typename T = Real>
class OpenBCSolver
{
public:
    using MF = typename R2C<T>::MF;
    using cMF = typename R2C<T>::cMF;

    /**
     * \brief Build a solver over \p domain using the FFT Info settings in \p info.
     *
     * \param domain Physical domain (cell or nodal indexing).
     * \param info   Optional FFT configuration. Batched solves are not supported.
     */
    explicit OpenBCSolver (Box const& domain, Info const& info = Info{});

    /**
     * \brief Populate the spectral Green's function used by subsequent solves.
     *
     * \tparam F Callable type satisfying `T(int i, int j, int k)`.
     * \param greens_function Functor that evaluates the physical-space Green's function.
     */
    template <class F>
    void setGreensFunction (F const& greens_function);

    /**
     * \brief Populate the spectral Green's function from a real-space Green's
     *        function given as a MultiFab.
     *
     * \p G must cover the lags 0..PaddedLength()-1 in each direction, that is
     * the Box from Domain().smallEnd() to Domain().smallEnd()+PaddedLength()-1;
     * its index type does not matter, only the indices do. Any BoxArray and
     * DistributionMapping are accepted. As with the functor version,
     * the Green's function is assumed even in every direction. Not supported
     * in twod_mode.
     *
     * \param G Green's function values at lags 0..PaddedLength()-1 in each direction.
     */
    void setGreensFunction (MF const& G);

    /**
     * \brief Populate the Green's function from its real, even spectral
     *        representation on an oversampled periodic grid.
     *
     * This serves Green's functions that are known in Fourier space rather
     * than in real space, such as the truncated free-space kernels of Vico,
     * Greengard & Ferrando (J. Comput. Phys. 323, 2016). Let M = 2*K be the
     * number of cells of a periodic grid, per direction, with the cell size of
     * the problem. \p ghat(i,j,k) must return the discrete Fourier transform,
     * on that grid, of the real-space Green's function g at wave-number
     * index (i,j,k), 0 <= i <= K[0] etc.; evenness in every direction supplies
     * the other indices. Then g = IDFT(ghat) is computed (as a DCT-I of the
     * (K+1)^d unique values) and used at lags 0..PaddedLength()-1, as in
     * setGreensFunction(MF const&).
     *
     * The caller chooses K. It must be at least PaddedLength()-1 in every
     * direction, and large enough that the periodic images of g do not reach
     * the lags used, which extend to N-1+ng for a domain of N points when
     * solve() fills ng ghost cells of phi. For a kernel truncated at radius
     * L, which must exceed the largest source-to-target distance including
     * ghost cells, this means 2*K[d] > N[d]-1+ng[d] + L/dx[d], plus a cell
     * or two of margin for the ringing of the band-limited kernel beyond L.
     * Fast FFT lengths for K are recommended (FFT::nextFastLen).
     *
     * The DCT work arrays, about 2*(K+1)^d reals on CPU and about twice that
     * on GPU including the transform scratch, exist only during this call.
     * Neither twod_mode nor a zero half_period in any direction is supported.
     *
     * \tparam F Callable type satisfying `T(int i, int j, int k)`.
     * \param ghat        Functor that evaluates the spectral Green's function.
     * \param half_period K, half the period of the oversampled grid, per direction.
     */
    template <class F>
    void setGreensFunctionFromSpectral (F const& ghat, IntVect const& half_period);

    /**
     * \brief Solve for \p phi given right-hand side \p rho.
     *
     * \param phi Solution MultiFab (overwritten in place).
     * \param rho Right-hand-side MultiFab.
     */
    void solve (MF& phi, MF const& rho);

    /**
     * \brief Access the physical domain this solver was built for.
     *
     * \return Original (un-extended) problem domain.
     */
    [[nodiscard]] Box const& Domain () const { return m_domain; }

    /**
     * \brief Access the one-sided padded length used to build the internal FFT domain.
     *
     * \return Padded length before the open-boundary convolution domain is doubled.
     */
    [[nodiscard]] IntVect const& PaddedLength () const { return m_padded_length; }

    /**
     * \brief Access the spectral Green's function held by this solver.
     *
     * This lets a caller keep its own set of Green's functions and move them into and
     * out of the solver, for instance with std::swap, instead of recomputing them with
     * setGreensFunction(). Since solve() only reads the spectral Green's function, one
     * moved in this way needs no further preparation.
     *
     * \note In 2D, this is an alias of the Green's function R2C's spectral data.
     *       Callers may write through it, but must not swap or otherwise re-seat it,
     *       because a later setGreensFunction() writes into the aliased array.
     *
     * \return The spectral Green's function.
     */
    [[nodiscard]] cMF& greensFunctionFFT () noexcept { return m_G_fft; }

    /**
     * \brief Access the spectral Green's function held by this solver.
     *
     * \return The spectral Green's function.
     */
    [[nodiscard]] cMF const& greensFunctionFFT () const noexcept { return m_G_fft; }

private:
    static IntVect make_padded_length (Box const& domain, Info const& info);
    static Box make_grown_domain (Box const& domain, IntVect const& padded_len,
                                  Info const& info);
    // A read-only view of src, with its boxes retyped to typ (keeping their
    // corners, unlike amrex::convert) and shifted by offset. No data is
    // copied; the const_cast is needed because a FAB cannot be built on a
    // pointer to const.
    static MF make_view (MF const& src, IndexType typ, IntVect const& offset);

    Box m_domain;
    Info m_info;
    IntVect m_padded_length;
    R2C<T> m_r2c;
    cMF m_G_fft;
    std::unique_ptr<R2C<T>> m_r2c_green;
};

template <typename T>
IntVect OpenBCSolver<T>::make_padded_length (Box const& domain, Info const& info)
{
    IntVect len = domain.length();
    int ndims = AMREX_SPACEDIM;
#if (AMREX_SPACEDIM == 3)
    if (info.twod_mode) { ndims = 2; }
#else
    amrex::ignore_unused(info);
#endif
    if (info.openbc_padding) {
        for (int idim = 0; idim < ndims; ++idim) {
            len[idim] = FFT::nextFastLen(len[idim], info.openbc_padding_nfactors);
        }
    }
    return len;
}

template <typename T>
Box OpenBCSolver<T>::make_grown_domain (Box const& domain, IntVect const& padded_len,
                                        Info const& info)
{
    IntVect len = padded_len;
    int ndims = AMREX_SPACEDIM;
#if (AMREX_SPACEDIM == 3)
    if (info.twod_mode) { ndims = 2; }
#else
    amrex::ignore_unused(info);
#endif
    for (int idim = 0; idim < ndims; ++idim) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
            len[idim] <= std::numeric_limits<int>::max()/2,
            "FFT::OpenBCSolver: padded domain length exceeds int range");
        len[idim] *= 2;
    }
    return Box(domain.smallEnd(), domain.smallEnd()+len-IntVect(1), domain.ixType());
}

template <typename T>
OpenBCSolver<T>::OpenBCSolver (Box const& domain, Info const& info)
    : m_domain(domain),
      m_info(info),
      m_padded_length(OpenBCSolver<T>::make_padded_length(domain, info)),
      m_r2c(OpenBCSolver<T>::make_grown_domain(domain, m_padded_length, info),
            m_info.setDomainStrategy(FFT::DomainStrategy::slab))
{
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(m_info.batch_size == 1,
        "FFT::OpenBCSolver does not support FFT::Info::batch_size > 1");

#if (AMREX_SPACEDIM == 3)
    if (m_info.twod_mode) {
        auto gdom = make_grown_domain(domain, m_padded_length, m_info);
        gdom.enclosedCells(2);
        gdom.setSmall(2, 0);
        int nprocs = std::min({ParallelContext::NProcsSub(),
                               m_info.nprocs,
                               m_domain.length(2)});
        gdom.setBig(2, nprocs-1);
        m_r2c_green = std::make_unique<R2C<T>>(gdom,m_info);
        auto [sd, ord] = m_r2c_green->getSpectralData();
        m_G_fft = cMF(*sd, amrex::make_alias, 0, 1);
    } else
#endif
    {
        amrex::ignore_unused(m_r2c_green);
        auto [sd, ord] = m_r2c.getSpectralData();
        amrex::ignore_unused(ord);
        m_G_fft.define(sd->boxArray(), sd->DistributionMap(), 1, 0);
    }
}

template <typename T>
template <class F>
void OpenBCSolver<T>::setGreensFunction (F const& greens_function)
{
    BL_PROFILE("OpenBCSolver::setGreensFunction");

    auto* infab = m_info.twod_mode ? detail::get_fab(m_r2c_green->m_rx)
        :                             detail::get_fab(m_r2c.m_rx);
    auto const& lo = m_domain.smallEnd();
    auto const& lo3 = lo.dim3();
    auto const len3d = m_padded_length.dim3();
    GpuArray<int,3> len{len3d.x, len3d.y, len3d.z};
    if (infab) {
        auto const& a = infab->array();
        auto box = infab->box();
        GpuArray<int,3> nimages{1,1,1};
        int ndims = m_info.twod_mode ? AMREX_SPACEDIM-1 : AMREX_SPACEDIM;
        for (int idim = 0; idim < ndims; ++idim) {
            if (box.smallEnd(idim) == lo[idim] && box.length(idim) == 2*len[idim]) {
                box.growHi(idim, -len[idim]+1); // +1 to include the middle plane
                nimages[idim] = 2;
            }
        }
        AMREX_ASSERT(nimages[0] == 2);
        box.shift(-lo);
        amrex::ParallelForOMP(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            T G;
            if (i == len[0] || j == len[1] || k == len[2]) {
                G = 0;
            } else {
                auto ii = i;
                auto jj = (j > len[1]) ? 2*len[1]-j : j;
                auto kk = (k > len[2]) ? 2*len[2]-k : k;
                G = greens_function(ii+lo3.x,jj+lo3.y,kk+lo3.z);
            }
            for (int koff = 0; koff < nimages[2]; ++koff) {
                int k2 = (koff == 0) ?  k : 2*len[2]-k;
                if ((k2 == 2*len[2]) || (koff == 1 && k == len[2])) {
                    continue;
                }
                for (int joff = 0; joff < nimages[1]; ++joff) {
                    int j2 = (joff == 0) ?  j : 2*len[1]-j;
                    if ((j2 == 2*len[1]) || (joff == 1 && j == len[1])) {
                        continue;
                    }
                    for (int ioff = 0; ioff < nimages[0]; ++ioff) {
                        int i2 = (ioff == 0) ?  i : 2*len[0]-i;
                        if ((i2 == 2*len[0]) || (ioff == 1 && i == len[0])) {
                            continue;
                        }
                        a(i2+lo3.x,j2+lo3.y,k2+lo3.z) = G;
                    }
                }
            }
        });
    }

    if (m_info.twod_mode) {
        m_r2c_green->forward(m_r2c_green->m_rx);
    } else {
        m_r2c.forward(m_r2c.m_rx);
    }

    if (!m_info.twod_mode) {
        auto [sd, ord] = m_r2c.getSpectralData();
        amrex::ignore_unused(ord);
        auto const* srcfab = detail::get_fab(*sd);
        if (srcfab) {
            auto* dstfab = detail::get_fab(m_G_fft);
            if (dstfab) {
                Gpu::dtod_memcpy_async(dstfab->dataPtr(), srcfab->dataPtr(), dstfab->nBytes());
            } else {
                amrex::Abort("FFT::OpenBCSolver: how did this happen");
            }
        }

        m_r2c.prepare_openbc();
    }
}

template <typename T>
auto OpenBCSolver<T>::make_view (MF const& src, IndexType typ, IntVect const& offset) -> MF
{
    using FAB = typename MF::FABType::value_type;
    auto transform = [&] (Box const& b) {
        return Box(b.smallEnd()+offset, b.bigEnd()+offset, typ);
    };
    BoxList bl(typ);
    for (int ibox = 0; ibox < src.boxArray().size(); ++ibox) {
        bl.push_back(transform(src.boxArray()[ibox]));
    }
    MF view(BoxArray(std::move(bl)), src.DistributionMap(), 1, src.nGrowVect(),
            MFInfo().SetAlloc(false));
    for (MFIter mfi(src); mfi.isValid(); ++mfi) {
        auto const& fab = src[mfi];
        view.setFab(mfi, FAB(transform(fab.box()), 1, const_cast<T*>(fab.dataPtr())));
    }
    return view;
}

template <typename T>
void OpenBCSolver<T>::setGreensFunction (MF const& G)
{
    BL_PROFILE("OpenBCSolver::setGreensFunction(MF)");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_info.twod_mode,
        "FFT::OpenBCSolver::setGreensFunction(MF): twod_mode is not supported");

    auto const& lo = m_domain.smallEnd();
    // lags 0..len-1; the plane at len is set to zero by the functor version
    Box const octant(lo, lo+m_padded_length-IntVect(1), m_domain.ixType());

    // Only the indices matter, not the index type.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(
        G.boxArray().contains(Box(octant.smallEnd(), octant.bigEnd(), G.ixType())),
        "FFT::OpenBCSolver::setGreensFunction(MF): G must cover the octant");

    // If G's index type differs from the domain's, view its data through
    // boxes of the domain's type.
    MF Gview;
    MF const* Gp = &G;
    if (G.ixType() != octant.ixType()) {
        Gview = make_view(G, octant.ixType(), IntVect(0));
        Gp = &Gview;
    }

    // Each process holds a part of the doubled domain, in which lag l and lag
    // 2*len-l share the value at lag l. Gather, per process, the hull of the
    // octant indices its part mirrors onto. m_rx has one box per process.
    auto const& rx = m_r2c.m_rx;
    BoxList bl(octant.ixType());
    for (int ibox = 0; ibox < rx.size(); ++ibox) {
        Box const& b = rx.box(ibox);
        IntVect mlo, mhi;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            int const len = m_padded_length[idim];
            int const a = b.smallEnd(idim) - lo[idim];
            int const c = b.bigEnd(idim) - lo[idim];
            int const ma = (a > len) ? 2*len-a : a;
            int const mc = (c > len) ? 2*len-c : c;
            // clamp to len-1: the plane at len is never read
            mhi[idim] = lo[idim] + std::min(len-1, (a <= len && len <= c) ? len : std::max(ma, mc));
            mlo[idim] = std::min(lo[idim] + std::min(ma, mc), mhi[idim]);
        }
        bl.push_back(Box(mlo, mhi, octant.ixType()));
    }
    MF gl(BoxArray(std::move(bl)), rx.DistributionMap(), 1, 0);
    gl.ParallelCopy(*Gp, 0, 0, 1);

    int const myproc = ParallelContext::MyProcSub();
    if (myproc < gl.size()) {
        auto const& ga = gl.const_array(myproc);
        this->setGreensFunction([=] AMREX_GPU_DEVICE (int i, int j, int k) -> T
        {
            return ga(i,j,k);
        });
    } else {
        // no local data; participate in the collective transform
        this->setGreensFunction([=] AMREX_GPU_DEVICE (int, int, int) -> T
        {
            return T(0);
        });
    }
    Gpu::streamSynchronize(); // gl goes out of scope
}

template <typename T>
template <class F>
void OpenBCSolver<T>::setGreensFunctionFromSpectral (F const& ghat, IntVect const& half_period)
{
    BL_PROFILE("OpenBCSolver::setGreensFunctionFromSpectral");

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!m_info.twod_mode,
        "FFT::OpenBCSolver::setGreensFunctionFromSpectral: twod_mode is not supported");
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(half_period[idim]+1 >= m_padded_length[idim],
            "FFT::OpenBCSolver::setGreensFunctionFromSpectral: "
            "half_period must be >= PaddedLength()-1");
    }

    // The spectral Green's function is real and even in every direction, so
    // of the (2K)^d samples on the periodic grid only (K+1)^d are unique, and
    // the inverse DFT reduces to a DCT-I. The samples are written directly
    // into the R2X's last-stage array and only its backward pass is run.
    // The R2X is local: this is called rarely compared to solve(), and its
    // work arrays would otherwise be held for the life of the solver.
    Box const dct_domain(IntVect(0), half_period);
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> bc;
    bc.fill(std::make_pair(Boundary::even_node, Boundary::even_node));
    Info info{};
    info.setNumProcs(m_info.nprocs);
    R2X<T> r2x(dct_domain, bc, info);

    // With all boundaries even_node, the last stage of the R2X is the real
    // array of the last dimension.
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(r2x.m_r2x_sub == nullptr,
        "FFT::OpenBCSolver::setGreensFunctionFromSpectral: "
        "half_period must be > 0 in every direction");
    T const scale = r2x.scalingFactor(); // 1/(2K)^d, the inverse DFT normalization
    auto fill = [=] AMREX_GPU_DEVICE (int i, int j, int k, auto& sp)
    {
        sp = scale * ghat(i,j,k);
    };
    // post_forward_doit maps storage indices of the transposed array to
    // logical (x,y,z) before calling fill.
#if (AMREX_SPACEDIM == 3)
    AMREX_ALWAYS_ASSERT(!r2x.m_rz.empty() && r2x.m_cz.empty());
    r2x.template post_forward_doit<2>(detail::get_fab(r2x.m_rz), fill);
#elif (AMREX_SPACEDIM == 2)
    AMREX_ALWAYS_ASSERT(!r2x.m_ry.empty() && r2x.m_cy.empty());
    r2x.template post_forward_doit<1>(detail::get_fab(r2x.m_ry), fill);
#else
    AMREX_ALWAYS_ASSERT(r2x.m_cx.empty());
    r2x.template post_forward_doit<0>(detail::get_fab(r2x.m_rx), fill);
#endif
    r2x.backward();

    // The R2X's real-space array now holds the Green's function at lags
    // 0..K from index 0; the MF version expects lags from Domain().smallEnd().
    this->setGreensFunction(make_view(r2x.m_rx, m_domain.ixType(), m_domain.smallEnd()));
}

template <typename T>
void OpenBCSolver<T>::solve (MF& phi, MF const& rho)
{
    BL_PROFILE("OpenBCSolver::solve");

    auto& inmf = m_r2c.m_rx;
    inmf.setVal(T(0));
    inmf.ParallelCopy(rho, 0, 0, 1);

    m_r2c.m_openbc_half = !m_info.twod_mode;
    m_r2c.forward(inmf);
    m_r2c.m_openbc_half = false;

    auto scaling_factor = m_r2c.scalingFactor();

    auto const* gfab = detail::get_fab(m_G_fft);
    if (gfab) {
        auto [sd, ord] = m_r2c.getSpectralData();
        amrex::ignore_unused(ord);
        auto* rhofab = detail::get_fab(*sd);
        if (rhofab) {
            auto* pdst = rhofab->dataPtr();
            auto const* psrc = gfab->dataPtr();
            Box const& rhobox = rhofab->box();
#if (AMREX_SPACEDIM == 3)
            Long leng = gfab->box().numPts();
            if (m_info.twod_mode) {
                AMREX_ASSERT(gfab->box().length(2) == 1 &&
                             leng == (rhobox.length(0) * rhobox.length(1)));
            } else {
                AMREX_ASSERT(leng == rhobox.numPts());
            }
#endif
            amrex::ParallelForOMP(rhobox.numPts(), [=] AMREX_GPU_DEVICE (Long i)
            {
#if (AMREX_SPACEDIM == 3)
                Long isrc = i % leng;
#else
                Long isrc = i;
#endif
                pdst[i] *= psrc[isrc] * scaling_factor;
            });
        } else {
            amrex::Abort("FFT::OpenBCSolver::solve: how did this happen?");
        }
    }

    m_r2c.m_openbc_half = !m_info.twod_mode;
    m_r2c.backward_doit(phi, phi.nGrowVect());
    m_r2c.m_openbc_half = false;
}

}

#endif
