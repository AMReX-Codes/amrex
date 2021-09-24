#ifndef CNS_DIFFUSION_K_H_
#define CNS_DIFFUSION_K_H_

#include "CNS_index_macros.H"
#include "CNS_parm.H"
#include <AMReX_FArrayBox.H>
#include <AMReX_CONSTANTS.H>
#include <cmath>

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diffcoef (int i, int j, int k,
              amrex::Array4<amrex::Real const> const& q,
              amrex::Array4<amrex::Real> const& coefs,
              Parm const& parm) noexcept
{
    using amrex::Real;

    if (q(i,j,k,QTEMP) < 0.)
    {
//        std::cout << "Computing diff coeffs at " << IntVect(i,j,k) <<
//             " with  temp = " << q(i,j,k,UTEMP)  << std::endl;
        amrex::Abort("Negative temperature sent to cns_diffcoef_eb");
    }

    coefs(i,j,k,CETA) = parm.C_S * std::sqrt(q(i,j,k,QTEMP)) * q(i,j,k,QTEMP) / (q(i,j,k,QTEMP)+parm.T_S);
    coefs(i,j,k,CXI)  = Real(0.0);
    coefs(i,j,k,CLAM) = coefs(i,j,k,CETA)*parm.cp/parm.Pr;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_constcoef (int i, int j, int k,
              amrex::Array4<amrex::Real> const& coefs,
              Parm const& parm) noexcept
{
    coefs(i,j,k,CETA) = parm.const_visc_mu;
    coefs(i,j,k,CXI)  = parm.const_visc_ki;
    coefs(i,j,k,CLAM) = parm.const_lambda;
}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_x (int i, int j, int k,
               amrex::Array4<amrex::Real const> const& q,
               amrex::Array4<amrex::Real const> const& coeffs,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
               amrex::Array4<amrex::Real> const& fx) noexcept
{
    using amrex::Real;

    Real  dTdx = (q(i,j,k,QTEMP)-q(i-1,j,k,QTEMP))*dxinv[0];
    AMREX_D_TERM(Real  dudx = (q(i,j,k,QU)-q(i-1,j,k,QU))*dxinv[0];,
                 Real  dvdx = (q(i,j,k,QV)-q(i-1,j,k,QV))*dxinv[0];,
                 Real  dwdx = (q(i,j,k,QW)-q(i-1,j,k,QW))*dxinv[0];);
    Real  dudy = (q(i,j+1,k,QU)+q(i-1,j+1,k,QU)-q(i,j-1,k,QU)-q(i-1,j-1,k,QU))*(0.25*dxinv[1]);
    Real  dvdy = (q(i,j+1,k,QV)+q(i-1,j+1,k,QV)-q(i,j-1,k,QV)-q(i-1,j-1,k,QV))*(0.25*dxinv[1]);
#if (AMREX_SPACEDIM == 2)
    Real  divu = dudx + dvdy;
#else
    Real  dudz = (q(i,j,k+1,QU)+q(i-1,j,k+1,QU)-q(i,j,k-1,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[2]);
    Real  dwdz = (q(i,j,k+1,QW)+q(i-1,j,k+1,QW)-q(i,j,k-1,QW)-q(i-1,j,k-1,QW))*(0.250*dxinv[2]);
    Real  divu = dudx + dvdy + dwdz;
#endif

    Real  etaf = 0.5*(coeffs(i,j,k,CETA)+coeffs(i-1,j,k,CETA));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i-1,j,k,CXI));

    AMREX_D_TERM(Real  tauxx = etaf*(2.0*dudx-(2.0/3.0)*divu) + xif*divu;,
                 Real  tauxy = etaf*(dudy+dvdx);,
                 Real  tauxz = etaf*(dudz+dwdx););

    AMREX_D_TERM(fx(i,j,k,UMX)   += -tauxx;,
                 fx(i,j,k,UMY)   += -tauxy;,
                 fx(i,j,k,UMZ)   += -tauxz;);
    fx(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i-1,j,k,QU))*tauxx+
                              (q(i,j,k,QV)+q(i-1,j,k,QV))*tauxy+
#if (AMREX_SPACEDIM == 3)
                              (q(i,j,k,QW)+q(i-1,j,k,QW))*tauxz+
#endif
                              (coeffs(i,j,k,CLAM) +coeffs(i-1,j,k,CLAM))*dTdx );

}

AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_y (int i, int j, int k, amrex::Array4<amrex::Real const> const& q,
               amrex::Array4<amrex::Real const> const& coeffs,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
               amrex::Array4<amrex::Real> const& fy) noexcept
{
    using amrex::Real;

    Real  dTdy = (q(i,j,k,QTEMP)-q(i,j-1,k,QTEMP))*dxinv[1];
    AMREX_D_TERM(Real  dudy = (q(i,j,k,QU)-q(i,j-1,k,QU))*dxinv[1];,
                 Real  dvdy = (q(i,j,k,QV)-q(i,j-1,k,QV))*dxinv[1];,
                 Real  dwdy = (q(i,j,k,QW)-q(i,j-1,k,QW))*dxinv[1];);

    Real  dudx = (q(i+1,j,k,QU)+q(i+1,j-1,k,QU)-q(i-1,j,k,QU)-q(i-1,j-1,k,QU))*(0.25*dxinv[0]);
    Real  dvdx = (q(i+1,j,k,QV)+q(i+1,j-1,k,QV)-q(i-1,j,k,QV)-q(i-1,j-1,k,QV))*(0.25*dxinv[0]);
#if (AMREX_SPACEDIM == 2)
    Real  divu = dudx + dvdy;
#else
    Real  dvdz = (q(i,j,k+1,QV)+q(i,j-1,k+1,QV)-q(i,j,k-1,QV)-q(i,j-1,k-1,QV))*(0.25*dxinv[2]);
    Real  dwdz = (q(i,j,k+1,QW)+q(i,j-1,k+1,QW)-q(i,j,k-1,QW)-q(i,j-1,k-1,QW))*(0.25*dxinv[2]);
    Real  divu = dudx + dvdy + dwdz;
#endif

    Real  etaf = 0.5*(coeffs(i,j,k,CETA)+coeffs(i,j-1,k,CETA));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j-1,k,CXI));
    AMREX_D_TERM(Real  tauyy = etaf*(2.0*dvdy-(2.0/3.0)*divu) + xif*divu;,
                 Real  tauxy = etaf*(dudy+dvdx);,
                 Real  tauyz = etaf*(dwdy+dvdz););

    AMREX_D_TERM(fy(i,j,k,UMX)   += -tauxy;,
                 fy(i,j,k,UMY)   += -tauyy;,
                 fy(i,j,k,UMZ)   += -tauyz;);
    fy(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i,j-1,k,QU))*tauxy
                             +(q(i,j,k,QV)+q(i,j-1,k,QV))*tauyy
#if (AMREX_SPACEDIM == 3)
                             +(q(i,j,k,QW)+q(i,j-1,k,QW))*tauyz
#endif
                             +(coeffs(i,j,k,CLAM) +coeffs(i,j-1,k,CLAM))*dTdy );

}

#if (AMREX_SPACEDIM == 3)
AMREX_GPU_DEVICE AMREX_FORCE_INLINE
void
cns_diff_z (int i, int j, int k,
               amrex::Array4<amrex::Real const> const& q,
               amrex::Array4<amrex::Real const> const& coeffs,
               amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxinv,
               amrex::Array4<amrex::Real> const& fz) noexcept
{
    using amrex::Real;

    Real  dTdz = (q(i,j,k,QTEMP)-q(i,j,k-1,QTEMP))*dxinv[2];
    Real  dudz = (q(i,j,k,QU)-q(i,j,k-1,QU))*dxinv[2];
    Real  dvdz = (q(i,j,k,QV)-q(i,j,k-1,QV))*dxinv[2];
    Real  dwdz = (q(i,j,k,QW)-q(i,j,k-1,QW))*dxinv[2];
    Real  dudx = (q(i+1,j,k,QU)+q(i+1,j,k-1,QU)-q(i-1,j,k,QU)-q(i-1,j,k-1,QU))*(0.25*dxinv[0]);
    Real  dwdx = (q(i+1,j,k,QW)+q(i+1,j,k-1,QW)-q(i-1,j,k,QW)-q(i-1,j,k-1,QW))*(0.25*dxinv[0]);
    Real  dvdy = (q(i,j+1,k,QV)+q(i,j+1,k-1,QV)-q(i,j-1,k,QV)-q(i,j-1,k-1,QV))*(0.25*dxinv[1]);
    Real  dwdy = (q(i,j+1,k,QW)+q(i,j+1,k-1,QW)-q(i,j-1,k,QW)-q(i,j-1,k-1,QW))*(0.25*dxinv[1]);
    Real  divu = dudx + dvdy + dwdz;
    Real  etaf = 0.5*(coeffs(i,j,k,CETA)+coeffs(i,j,k-1,CETA));
    Real  xif  = 0.5*(coeffs(i,j,k,CXI)+coeffs(i,j,k-1,CXI));
    Real  tauxz = etaf*(dudz+dwdx);
    Real  tauyz = etaf*(dvdz+dwdy);
    Real  tauzz = etaf*(2.0*dwdz-(2.0/3.0)*divu) + xif*divu;

    fz(i,j,k,UMX)   += -tauxz;
    fz(i,j,k,UMY)   += -tauyz;
    fz(i,j,k,UMZ)   += -tauzz;
    fz(i,j,k,UEDEN) += -0.5*( (q(i,j,k,QU)+q(i,j,k-1,QU))*tauxz
                             +(q(i,j,k,QV)+q(i,j,k-1,QV))*tauyz
                             +(q(i,j,k,QW)+q(i,j,k-1,QW))*tauzz
                             +(coeffs(i,j,k,CLAM) +coeffs(i,j,k-1,CLAM))*dTdz );
}
#endif

#endif
