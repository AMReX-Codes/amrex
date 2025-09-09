#ifndef AMREX_FFT_OPENBC_SOLVER_H_
#define AMREX_FFT_OPENBC_SOLVER_H_

#include <AMReX_FFT_R2C.H>

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

private:
    static IntVect make_padded_length (Box const& domain, Info const& info);
    static Box make_grown_domain (Box const& domain, IntVect const& padded_len,
                                  Info const& info);

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
        amrex::ParallelForOMP(box, [=]<int WIDTH> AMREX_GPU_DEVICE (amrex::SIMDindex<WIDTH> i, int j, int k)
        {
            using SIMD_T = simd::stdx::fixed_size_simd<T, i.width>;
            using SIMD_int = simd::stdx::fixed_size_simd<int, i.width>; // simd::stdx::rebind_simd_t<int, SIMD_T>;

            SIMD_T G = 0;
            if (j != len[1] && k != len[2])
            {
                SIMD_int ii = simd::stdx::iota_v<SIMD_int> + i.index;
                auto jj = (j > len[1]) ? 2*len[1]-j : j;
                auto kk = (k > len[2]) ? 2*len[2]-k : k;

                auto const i_bound = ii == len[0];
                simd::stdx::where(simd::stdx::cvt(i_bound), G) = 0.0;
                simd::stdx::where(simd::stdx::cvt(!i_bound), G) = greens_function.template operator()<WIDTH>(ii+lo3.x,jj+lo3.y,kk+lo3.z);
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
                        for (int iw = i.index; iw < i.index+i.width; ++iw) {
                            int i2 = (ioff == 0) ?  iw : 2*len[0]-iw;
                            if ((i2 == 2*len[0]) || (ioff == 1 && iw == len[0])) {
                                continue;
                            }
                            // TODO: SIMD-assign N values
                            a(i2+lo3.x,j2+lo3.y,k2+lo3.z) = G[iw];
                        }
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
