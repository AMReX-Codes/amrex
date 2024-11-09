#ifndef AMREX_FFT_OPENBC_SOLVER_H_
#define AMREX_FFT_OPENBC_SOLVER_H_

#include <AMReX_FFT_R2C.H>

namespace amrex::FFT
{

template <typename T = Real>
class OpenBCSolver
{
public:
    using MF = typename R2C<T>::MF;
    using cMF = typename R2C<T>::cMF;

    explicit OpenBCSolver (Box const& domain, Info const& info = Info{});

    template <class F>
    void setGreensFunction (F const& greens_function);

    void solve (MF& phi, MF const& rho);

    [[nodiscard]] Box const& Domain () const { return m_domain; }

private:
    static Box make_grown_domain (Box const& domain, Info const& info);

    Box m_domain;
    Info m_info;
    R2C<T> m_r2c;
    cMF m_G_fft;
    std::unique_ptr<R2C<T>> m_r2c_green;
};

template <typename T>
Box OpenBCSolver<T>::make_grown_domain (Box const& domain, Info const& info)
{
    IntVect len = domain.length();
#if (AMREX_SPACEDIM == 3)
    if (info.batch_mode) { len[2] = 0; }
#else
    amrex::ignore_unused(info);
#endif
    return Box(domain.smallEnd(), domain.bigEnd()+len, domain.ixType());
}

template <typename T>
OpenBCSolver<T>::OpenBCSolver (Box const& domain, Info const& info)
    : m_domain(domain),
      m_info(info),
      m_r2c(OpenBCSolver<T>::make_grown_domain(domain,info), info)
{
#if (AMREX_SPACEDIM == 3)
    if (m_info.batch_mode) {
#if defined(AMREX_USE_SYCL)
        amrex::Abort("FFT::OpenBCSolver: batch mode not support because Intel does not support batched 2d fft.");
#endif
        auto gdom = make_grown_domain(domain,m_info);
        gdom.enclosedCells(2);
        gdom.setSmall(2, 0);
        int nprocs = std::min({ParallelContext::NProcsSub(),
                               m_info.nprocs,
                               m_domain.length(2)});
        gdom.setBig(2, nprocs-1);
        m_r2c_green = std::make_unique<R2C<T>>(gdom,info);
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
    auto* infab = m_info.batch_mode ? detail::get_fab(m_r2c_green->m_rx)
        :                             detail::get_fab(m_r2c.m_rx);
    auto const& lo = m_domain.smallEnd();
    auto const& lo3 = lo.dim3();
    auto const& len = m_domain.length3d();
    if (infab) {
        auto const& a = infab->array();
        auto box = infab->box();
        GpuArray<int,3> nimages{1,1,1};
        int ndims = m_info.batch_mode ? AMREX_SPACEDIM : AMREX_SPACEDIM-1;
        for (int idim = 0; idim < ndims; ++idim) {
            if (box.smallEnd(idim) == lo[idim] && box.length(idim) == 2*len[idim]) {
                box.growHi(idim, -len[idim]+1); // +1 to include the middle plane
                nimages[idim] = 2;
            }
        }
        AMREX_ASSERT(nimages[0] == 2);
        box.shift(-lo);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
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

    if (m_info.batch_mode) {
        m_r2c_green->forward(m_r2c_green->m_rx);
    } else {
        m_r2c.forward(m_r2c.m_rx);
    }

    if (!m_info.batch_mode) {
        auto [sd, ord] = m_r2c.getSpectralData();
        amrex::ignore_unused(ord);
        auto const* srcfab = detail::get_fab(*sd);
        if (srcfab) {
            auto* dstfab = detail::get_fab(m_G_fft);
            if (dstfab) {
#if defined(AMREX_USE_GPU)
                Gpu::dtod_memcpy_async
#else
                std::memcpy
#endif
                    (dstfab->dataPtr(), srcfab->dataPtr(), dstfab->nBytes());
            } else {
                amrex::Abort("FFT::OpenBCSolver: how did this happen");
            }
        }
    }
}

template <typename T>
void OpenBCSolver<T>::solve (MF& phi, MF const& rho)
{
    auto& inmf = m_r2c.m_rx;
    inmf.setVal(T(0));
    inmf.ParallelCopy(rho, 0, 0, 1);

    m_r2c.forward(inmf);

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
            if (m_info.batch_mode) {
                AMREX_ASSERT(gfab->box().length(2) == 1 &&
                             leng == (rhobox.length(0) * rhobox.length(1)));
            } else {
                AMREX_ASSERT(leng == rhobox.numPts());
            }
#endif
            amrex::ParallelFor(rhobox.numPts(), [=] AMREX_GPU_DEVICE (Long i)
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

    m_r2c.backward_doit(phi, phi.nGrowVect());
}

}

#endif
