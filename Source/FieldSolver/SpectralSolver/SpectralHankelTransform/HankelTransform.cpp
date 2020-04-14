/* Copyright 2019 David Grote
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#include "WarpX.H"

#include "HankelTransform.H"
#include "BesselRoots.H"

#include <blas.hh>
#include <lapack.hh>

using amrex::operator""_rt;

HankelTransform::HankelTransform (int const hankel_order,
                                  int const azimuthal_mode,
                                  int const nr,
                                  const amrex::Real rmax)
: m_nr(nr), m_nk(nr)
{

    // Check that azimuthal_mode has a valid value
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(hankel_order-1 <= azimuthal_mode && azimuthal_mode <= hankel_order+1,
                                     "azimuthal_mode must be either hankel_order-1, hankel_order or hankel_order+1");

    RealVector alphas;
    amrex::Vector<int> alpha_errors;

    GetBesselRoots(azimuthal_mode, m_nk, alphas, alpha_errors);
    // Add of check of alpha_errors, should be all zeros
    AMREX_ALWAYS_ASSERT(std::all_of(alpha_errors.begin(), alpha_errors.end(), [](int i) { return i == 0; }));

    // Calculate the spectral grid
    m_kr.resize(m_nk);
    for (int ik=0 ; ik < m_nk ; ik++) {
        m_kr[ik] = alphas[ik]/rmax;
    }

    // Calculate the spatial grid (Uniform grid with a half-cell offset)
    RealVector rmesh(m_nr);
    amrex::Real dr = rmax/m_nr;
    for (int ir=0 ; ir < m_nr ; ir++) {
        rmesh[ir] = dr*(ir + 0.5_rt);
    }

    // Calculate and store the inverse matrix invM
    // (imposed by the constraints on the DHT of Bessel modes)
    // NB: When compared with the FBPIC article, all the matrices here
    // are calculated in transposed form. This is done so as to use the
    // `dot` and `gemm` functions, in the `transform` method.
    int p_denom;
    if (hankel_order == azimuthal_mode) {
        p_denom = hankel_order + 1;
    } else {
        p_denom = hankel_order;
    }

    RealVector denom(m_nk);
    for (int ik=0 ; ik < m_nk ; ik++) {
        const amrex::Real jna = jn(p_denom, alphas[ik]);
        denom[ik] = MathConst::pi*rmax*rmax*jna*jna;
    }

    RealVector num(m_nk*m_nr);
    for (int ir=0 ; ir < m_nr ; ir++) {
        for (int ik=0 ; ik < m_nk ; ik++) {
            int const ii = ik + ir*m_nk;
            num[ii] = jn(hankel_order, rmesh[ir]*m_kr[ik]);
        }
    }

    // Get the inverse matrix
    invM.resize(m_nk*m_nr);
    if (azimuthal_mode > 0) {
        for (int ir=0 ; ir < m_nr ; ir++) {
            for (int ik=1 ; ik < m_nk ; ik++) {
                int const ii = ik + ir*m_nk;
                invM[ii] = num[ii]/denom[ik];
            }
        }
        // ik = 0
        // In this case, the functions are represented by Bessel functions
        // *and* an additional mode (below) which satisfies the same
        // algebric relations for curl/div/grad as the regular Bessel modes,
        // with the value kperp=0.
        // The normalization of this mode is arbitrary, and is chosen
        // so that the condition number of invM is close to 1
        if (hankel_order == azimuthal_mode-1) {
            for (int ir=0 ; ir < m_nr ; ir++) {
                int const ii = ir*m_nk;
                invM[ii] = std::pow(rmesh[ir], (azimuthal_mode-1))/(MathConst::pi*std::pow(rmax, (azimuthal_mode+1)));
            }
        } else {
            for (int ir=0 ; ir < m_nr ; ir++) {
                int const ii = ir*m_nk;
                invM[ii] = 0.;
            }
        }
    } else {
        for (int ir=0 ; ir < m_nr ; ir++) {
            for (int ik=0 ; ik < m_nk ; ik++) {
                int const ii = ik + ir*m_nk;
                invM[ii] = num[ii]/denom[ik];
            }
        }
    }

    // Calculate the matrix M by inverting invM
    if (azimuthal_mode !=0 && hankel_order != azimuthal_mode-1) {
        // In this case, invM is singular, thus we calculate the pseudo-inverse.
        // The Moore-Penrose psuedo-inverse is calculated using the SVD method.

        M.resize(m_nk*m_nr, 0.);
        RealVector invMcopy(invM);
        RealVector sdiag(m_nk-1, 0.);
        RealVector u((m_nk-1)*(m_nk-1), 0.);
        RealVector vt((m_nr)*(m_nr), 0.);
        RealVector sp((m_nr)*(m_nk-1), 0.);
        RealVector temp((m_nr)*(m_nk-1), 0.);

        // Calculate the singlular-value-decomposition of invM (leaving out the first row).
        // invM = u*sdiag*vt
        // Note that invMcopy.dataPtr()+1 is passed in so that the first ik row is skipped
        // A copy is passed in since the matrix is destroyed
        lapack::gesvd(lapack::Job::AllVec, lapack::Job::AllVec,
                      m_nk-1, m_nr, invMcopy.dataPtr()+1, m_nk,
                      sdiag.dataPtr(), u.dataPtr(), m_nk-1,
                      vt.dataPtr(), m_nr);

        // Calculate the pseudo-inverse of sdiag, which is trivial since it
        // only has values of the diagonal.
        for (int i=0 ; i < m_nk-1 ; i++) {
            if (sdiag[i] != 0.) {
                int const j = i + i*m_nk;
                sp[j] = 1._rt/sdiag[i];
            }
        }

        // M can be calculated from the SVD
        // M = v*sp*u_transpose

        // Do the second multiply, temp = sp*u_transpose
        // temp(1:n,1:m) = matmul(transpose(vt(1:n,1:n)), matmul(sp(1:n,1:m), transpose(u(1:m,1:m))))
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans,
                   m_nr, m_nk-1, m_nk-1, 1._rt,
                   sp.dataPtr(), m_nr,
                   u.dataPtr(), m_nk-1, 0._rt,
                   temp.dataPtr(), m_nr);

        // Do the first mutiply, M = vt*temp
        // Note that M.dataPtr()+m_nr is passed in so that the first ir column is skipped
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
                   m_nr, m_nk-1, m_nr, 1.,
                   vt.dataPtr(), m_nr,
                   temp.dataPtr(), m_nr, 0._rt,
                   M.dataPtr()+m_nr, m_nr);

    } else {
        // In this case, invM is invertible; calculate the inverse.
        // getrf calculates the LU decomposition
        // getri calculates the inverse from the LU decomposition

        M = invM;
        amrex::Vector<int64_t> ipiv(m_nr);
        lapack::getrf(m_nk, m_nr, M.dataPtr(), m_nk, ipiv.dataPtr());
        lapack::getri(m_nr, M.dataPtr(), m_nr, ipiv.dataPtr());

    }

}

void
HankelTransform::HankelForwardTransform (amrex::FArrayBox const& F, int const F_icomp,
                                         amrex::FArrayBox      & G, int const G_icomp)
{
    amrex::Box const& F_box = F.box();
    amrex::Box const& G_box = G.box();

    int const nrF = F_box.length(0);
    int const nz = F_box.length(1);
    int const ngr = G_box.smallEnd(0) - F_box.smallEnd(0);

    AMREX_ALWAYS_ASSERT(m_nr == G_box.length(0));
    AMREX_ALWAYS_ASSERT(nz == G_box.length(1));
    AMREX_ALWAYS_ASSERT(ngr >= 0);
    AMREX_ALWAYS_ASSERT(F_box.bigEnd(0)+1 >= m_nr);

#ifndef AMREX_USE_GPU
    // On CPU, the blas::gemm is significantly faster

    // Note that M is flagged to be transposed since it has dimensions (m_nr, m_nk)
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
               m_nk, nz, m_nr, 1._rt,
               M.dataPtr(), m_nk,
               F.dataPtr(F_icomp)+ngr, nrF, 0._rt,
               G.dataPtr(G_icomp), m_nk);

#else
    // On GPU, the explicit loop is significantly faster
    // It is not clear if the GPU gemm wasn't build properly, it is cycling data out and back
    // in to the device, or if it is because gemm is launching its own threads.

    amrex::Real const * M_arr = M.dataPtr();
    amrex::Array4<const amrex::Real> const & F_arr = F.array();
    amrex::Array4<      amrex::Real> const & G_arr = G.array();

    int const nr = m_nr;

    ParallelFor(G_box,
    [=] AMREX_GPU_DEVICE(int ik, int iz, int inotused) noexcept {
        G_arr(ik,iz,G_icomp) = 0.;
        for (int ir=0 ; ir < nr ; ir++) {
            int const ii = ir + ik*nr;
            G_arr(ik,iz,G_icomp) += M_arr[ii]*F_arr(ir,iz,F_icomp);
        }
    });

#endif

}

void
HankelTransform::HankelInverseTransform (amrex::FArrayBox const& G, int const G_icomp,
                                         amrex::FArrayBox      & F, int const F_icomp)
{
    amrex::Box const& G_box = G.box();
    amrex::Box const& F_box = F.box();

    int const nrF = F_box.length(0);
    int const nz = F_box.length(1);
    int const ngr = G_box.smallEnd(0) - F_box.smallEnd(0);

    AMREX_ALWAYS_ASSERT(m_nr == G_box.length(0));
    AMREX_ALWAYS_ASSERT(nz == G_box.length(1));
    AMREX_ALWAYS_ASSERT(ngr >= 0);
    AMREX_ALWAYS_ASSERT(F_box.bigEnd(0)+1 >= m_nr);

#ifndef AMREX_USE_GPU
    // On CPU, the blas::gemm is significantly faster

    // Note that invM is flagged to be transposed since it has dimensions (m_nk, m_nr)
    blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
               m_nr, nz, m_nk, 1._rt,
               invM.dataPtr(), m_nr,
               G.dataPtr(G_icomp), m_nk, 0._rt,
               F.dataPtr(F_icomp)+ngr, nrF);

#else
    // On GPU, the explicit loop is significantly faster
    // It is not clear if the GPU gemm wasn't build properly, it is cycling data out and back
    // in to the device, or if it is because gemm is launching its own threads.

    amrex::Real const * invM_arr = invM.dataPtr();
    amrex::Array4<const amrex::Real> const & G_arr = G.array();
    amrex::Array4<      amrex::Real> const & F_arr = F.array();

    int const nk = m_nk;

    ParallelFor(G_box,
    [=] AMREX_GPU_DEVICE(int ir, int iz, int inotused) noexcept {
        F_arr(ir,iz,F_icomp) = 0.;
        for (int ik=0 ; ik < nk ; ik++) {
            int const ii = ik + ir*nk;
            F_arr(ir,iz,F_icomp) += invM_arr[ii]*G_arr(ik,iz,G_icomp);
        }
    });

#endif

}
