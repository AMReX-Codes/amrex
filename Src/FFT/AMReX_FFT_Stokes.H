#ifndef AMREX_FFT_STOKES_H_
#define AMREX_FFT_STOKES_H_

#include <AMReX_FFT.H>
#include <AMReX_Geometry.H>

namespace amrex::FFT
{

class Stokes
{
public:

    static_assert(AMREX_SPACEDIM >= 2, "FFT::Stokes requires 2D or 3D");

    /**
     * \brief Construct a Stokes solver with explicit boundary types.
     *
     * Only periodic BCs are supported.
     *
     * \param geom Geometry describing the domain and metric information.
     * \param bc   Pair of boundary descriptors (low/high) per coordinate direction.
     */
    Stokes (Geometry const& geom,
             Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> const& bc)
        : m_domain_lo(geom.Domain().smallEnd()),
          m_geom(detail::shift_geom(geom)),
          m_bc(bc)
    {
        bool all_periodic = true;
        for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
            all_periodic = all_periodic
                && (bc[idim].first == Boundary::periodic)
                && (bc[idim].second == Boundary::periodic);
        }
        if (all_periodic) {
            m_r2c = std::make_unique<R2C<Real>>(m_geom.Domain());
        } else {
            amrex::Abort("FFT::Stokes: only supports periodic BC");
        }
    }

    /**
     * \brief Construct a purely periodic Stokes solver.
     *
     * \param geom Periodic geometry (all directions must be periodic).
     */
    explicit Stokes (Geometry const& geom)
        : m_domain_lo(geom.Domain().smallEnd()),
          m_geom(detail::shift_geom(geom)),
          m_bc{AMREX_D_DECL(std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic),
                            std::make_pair(Boundary::periodic,Boundary::periodic))}
    {
        if (m_geom.isAllPeriodic()) {
            m_r2c = std::make_unique<R2C<Real>>(m_geom.Domain());
        } else {
            amrex::Abort("FFT::Stokes: only supports periodic BC");
        }
    }

    /**
     * \brief Solve the generalized Stokes problem in spectral space.
     *
     * Solves \f$(\alpha - \eta \nabla^2) \mathbf{u} + \nabla p = \mathbf{f}\f$
     * with the divergence-free constraint enforced in Fourier space.
     * This solver is for 2D and 3D only.
     *
     * This solver uses a MAC staggered grid. \p U and \p rhsx must be
     * at x-face centers, \p V and \p rhsy at y-face centers,
     * \p W and \p rhsz at z-face centers, and \p p cell-centered.
     *
     * \param U     x-component of velocity (output, at x-face centers).
     * \param V     y-component of velocity (output, at y-face centers).
     * \param W     z-component of velocity (output, at z-face centers, 3D only).
     * \param p     Pressure (output, cell-centered).
     * \param rhsx  x-component of the right-hand side (at x-face centers).
     * \param rhsy  y-component of the right-hand side (at y-face centers).
     * \param rhsz  z-component of the right-hand side (at z-face centers, 3D only).
     * \param geom  Geometry describing the domain.
     * \param alpha Scalar coefficient in the Helmholtz operator.
     * \param eta   Scalar coefficient in the Helmholtz operator.
     */
    void solve (AMREX_D_DECL(MultiFab& U, MultiFab& V, MultiFab& W), MultiFab& p,
                AMREX_D_DECL(MultiFab const& rhsx, MultiFab const& rhsy, MultiFab const& rhsz),
                Geometry const& geom, Real alpha, Real eta);

private:
    IntVect m_domain_lo;
    Geometry m_geom;
    Array<std::pair<Boundary,Boundary>,AMREX_SPACEDIM> m_bc;
    std::unique_ptr<R2C<Real>> m_r2c;
};


void Stokes::solve (MultiFab& U,
                    MultiFab& V,
#if (BL_SPACEDIM == 3)
                    MultiFab& W,
#endif
                    MultiFab& p,
                    MultiFab const& rhsx,
                    MultiFab const& rhsy,
#if (BL_SPACEDIM == 3)
                    MultiFab const& rhsz,
#endif
                    Geometry const& geom,
                    Real alpha, Real eta)
{

            auto const& dx = m_geom.CellSizeArray();
            auto scaling = Real(1) / Real(m_geom.Domain().d_numPts());

            R2C<Real> r2c(m_geom.Domain());
            auto const& [cba, cdm] = r2c.getSpectralDataLayout();

            cMultiFab phat(cba, cdm, 1, 0);

            cMultiFab uhat(cba,cdm,1,0);
            cMultiFab vhat(cba,cdm,1,0);
#if (BL_SPACEDIM ==3)
            cMultiFab what(cba,cdm,1,0);
#endif

            cMultiFab rxhat(cba,cdm,1,0);
            r2c.forward(rhsx, rxhat);
            amrex::Scale(rxhat, scaling, 0, 1, 0);

            cMultiFab ryhat(cba,cdm,1,0);
            r2c.forward(rhsy, ryhat);
            amrex::Scale(ryhat, scaling, 0, 1, 0);

#if (BL_SPACEDIM ==3)
            cMultiFab rzhat(cba,cdm,1,0);
            r2c.forward(rhsz, rzhat);
            amrex::Scale(rzhat, scaling, 0, 1, 0);
#endif

            using Complex = GpuComplex<Real>;
            int const nx = m_geom.Domain().length(0);
#if (AMREX_SPACEDIM >= 2)
            int const ny = m_geom.Domain().length(1);
#endif
#if (AMREX_SPACEDIM == 3)
            int const nz = m_geom.Domain().length(2);
#endif
            for (MFIter mfi(phat); mfi.isValid(); ++mfi) {
                auto const& pb = phat[mfi].box();
                auto const& p = phat[mfi].array();
                auto const& u = uhat[mfi].array();
                auto const& v = vhat[mfi].array();
#if (BL_SPACEDIM ==3)
                auto const& w = what[mfi].array();
#endif
                auto const& rx = rxhat[mfi].array();
                auto const& ry = ryhat[mfi].array();
#if (BL_SPACEDIM ==3)
                auto const& rz = rzhat[mfi].array();
#endif
                ParallelFor(pb, [=, this] AMREX_GPU_DEVICE (int i, int j, int k) noexcept {

                   AMREX_D_TERM(int ik = (i <= nx/2) ? i : i - nx;,
                                int jk = (j <= ny/2) ? j : j - ny;,
                                int kk = (k <= nz/2) ? k : k - nz);

                   GpuArray<Real,AMREX_SPACEDIM> kwave
                         {AMREX_D_DECL(Real(2)*ik*Math::pi<Real>()/Real(m_geom.Domain().length(0)),
                                       Real(2)*jk*Math::pi<Real>()/Real(m_geom.Domain().length(1)),
                                       Real(2)*kk*Math::pi<Real>()/Real(m_geom.Domain().length(2)))};

                   Real delsqk =
                         AMREX_D_TERM(Real(2)*(std::cos(kwave[0])-Real(1))/(dx[0]*dx[0]),
                                     +(Real(2)*(std::cos(kwave[1])-Real(1))/(dx[1]*dx[1])),
                                     +(Real(2)*(std::cos(kwave[2])-Real(1))/(dx[2]*dx[2])));

                   amrex::GpuArray<amrex::GpuComplex<amrex::Real>, AMREX_SPACEDIM> delkp
                         {AMREX_D_DECL(
                               amrex::GpuComplex<amrex::Real>((std::cos(kwave[0])-Real(1))/dx[0],
                                                               std::sin(kwave[0])/dx[0]),
                               amrex::GpuComplex<amrex::Real>((std::cos(kwave[1])-Real(1))/dx[1],
                                                               std::sin(kwave[1])/dx[1]),
                               amrex::GpuComplex<amrex::Real>((std::cos(kwave[2])-Real(1))/dx[2],
                                                               std::sin(kwave[2])/dx[2]))};

                   amrex::GpuArray<amrex::GpuComplex<amrex::Real>, AMREX_SPACEDIM> delkm
                         {AMREX_D_DECL(
                               amrex::GpuComplex<amrex::Real>((Real(1)-std::cos(kwave[0]))/dx[0],
                                                                  std::sin(kwave[0])/dx[0]),
                               amrex::GpuComplex<amrex::Real>((Real(1)-std::cos(kwave[1]))/dx[1],
                                                                  std::sin(kwave[1])/dx[1]),
                               amrex::GpuComplex<amrex::Real>((Real(1)-std::cos(kwave[2]))/dx[2],
                                                                  std::sin(kwave[2])/dx[2]))};

                   Complex rhsdotdp =
                         AMREX_D_TERM(rx(i,j,k)*delkp[0],
                                     +ry(i,j,k)*delkp[1],
                                     +rz(i,j,k)*delkp[2]);

                   if(std::abs(delsqk) > Real(1.e-10)){
                       p(i,j,k) = rhsdotdp / delsqk;
                   } else {
                       p(i,j,k) = Real(0.);
                   }

                   Real diffop = alpha - eta*delsqk;
                   if(diffop > Real(1.e-10)) {
                        u(i,j,k) = (rx(i,j,k) - p(i,j,k)*delkm[0])/diffop;
                        v(i,j,k) = (ry(i,j,k) - p(i,j,k)*delkm[1])/diffop;
#if (BL_SPACEDIM ==3)
                        w(i,j,k) = (rz(i,j,k) - p(i,j,k)*delkm[2])/diffop;
#endif
                   } else {
                        u(i,j,k) = Real(0.);
                        v(i,j,k) = Real(0.);
#if (BL_SPACEDIM ==3)
                        w(i,j,k) = Real(0.);
#endif
                   }

                 });
            }


                r2c.backward(phat, p);
                r2c.backward(uhat, U);
                r2c.backward(vhat, V);
#if (BL_SPACEDIM ==3)
                r2c.backward(what, W);
#endif

}

} // namespace amrex::FFT

#endif
