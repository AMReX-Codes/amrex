#ifndef AMREX_FFT_OPENBC_SOLVER_H_
#define AMREX_FFT_OPENBC_SOLVER_H_

#include <AMReX_FFT_R2C.H>

#include <AMReX_VisMF.H>

namespace amrex::FFT
{

template <typename T = Real>
class OpenBCSolver
{
public:
    using MF = typename R2C<T>::MF;
    using cMF = typename R2C<T>::cMF;

    explicit OpenBCSolver (Box const& domain);

    template <class F>
    void setGreensFunction (F const& greens_function);

    void solve (MF& phi, MF const& rho);

    [[nodiscard]] Box const& Domain () const { return m_domain; }

private:
    Box m_domain;
    R2C<T> m_r2c;
    cMF m_G_fft;
};

template <typename T>
OpenBCSolver<T>::OpenBCSolver (Box const& domain)
    : m_domain(domain),
      m_r2c(Box(domain.smallEnd(), domain.bigEnd()+domain.length(), domain.ixType()))
{
    auto [sd, ord] = m_r2c.getSpectralData();
    amrex::ignore_unused(ord);
    m_G_fft.define(sd->boxArray(), sd->DistributionMap(), 1, 0);
}

template <typename T>
template <class F>
void OpenBCSolver<T>::setGreensFunction (F const& greens_function)
{
    auto* infab = detail::get_fab(m_r2c.m_rx);
    auto const& lo = m_domain.smallEnd();
    auto const& lo3 = lo.dim3();
    auto const& len = m_domain.length3d();
    if (infab) {
        auto const& a = infab->array();
        auto box = infab->box();
        GpuArray<int,3> nimages{1,1,1};
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            if (box.smallEnd(idim) == lo[idim] && box.length(idim) == 2*len[idim]) {
                box.growHi(idim, -len[idim]+1); // +1 to include the middle plane
                nimages[idim] = 2;
            }
        }
        AMREX_ASSERT(nimages[0] == 2);
        box.shift(-lo);
        amrex::ParallelFor(box, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (i == len[0] || j == len[1] || k == len[2]) {
                a(i+lo3.x,j+lo3.y,k+lo3.z) = T(0);
            } else {
                auto ii = i;
                auto jj = (j > len[1]) ? 2*len[1]-j : j;
                auto kk = (k > len[2]) ? 2*len[2]-k : k;
                auto G = greens_function(ii+lo3.x,jj+lo3.y,kk+lo3.z);
                for (int koff = 0; koff < nimages[2]; ++koff) {
                    int k2 = (koff == 0) ?  k : 2*len[2]-k;
                    if (k2 == 2*len[2]) { continue; }
                    for (int joff = 0; joff < nimages[1]; ++joff) {
                        int j2 = (joff == 0) ?  j : 2*len[1]-j;
                        if (j2 == 2*len[1]) { continue; }
                        for (int ioff = 0; ioff < nimages[0]; ++ioff) {
                            int i2 = (ioff == 0) ?  i : 2*len[0]-i;
                            if (i2 == 2*len[0]) { continue; }
                            a(i2+lo3.x,j2+lo3.y,k2+lo3.z) = G;
                        }
                    }
                }
            }
        });
    }

    m_r2c.forward(m_r2c.m_rx);

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

template <typename T>
void OpenBCSolver<T>::solve (MF& phi, MF const& rho)
{
    auto& inmf = m_r2c.m_rx;
    inmf.setVal(T(0));
    inmf.ParallelCopy(rho, 0, 0, 1);

    m_r2c.forward(inmf);

    auto scaling_factor = T(1) / T(m_r2c.m_real_domain.numPts());

    auto const* gfab = detail::get_fab(m_G_fft);
    if (gfab) {
        auto [sd, ord] = m_r2c.getSpectralData();
        amrex::ignore_unused(ord);
        auto* rhofab = detail::get_fab(*sd);
        if (rhofab) {
            auto* pdst = rhofab->dataPtr();
            auto const* psrc = gfab->dataPtr();
            amrex::ParallelFor(rhofab->box().numPts(), [=] AMREX_GPU_DEVICE (Long i)
            {
                pdst[i] *= psrc[i] * scaling_factor;
            });
        } else {
            amrex::Abort("FFT::OpenBCSolver::solve: how did this happen?");
        }
    }

    m_r2c.backward_doit(phi, phi.nGrowVect());
}

}

#endif
