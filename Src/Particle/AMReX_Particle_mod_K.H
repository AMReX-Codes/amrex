#ifndef AMREX_PARTICLE_MOD_K_H_
#define AMREX_PARTICLE_MOD_K_H_

#include <AMReX_Particles.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

AMREX_GPU_HOST_DEVICE AMREX_INLINE
template <typename P>
void amrex_deposit_cic (P const& p, int nc, amrex::Array4<amrex::Real const> const& rho,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

    int i = std::floor(lx);
    int j = std::floor(ly);
    int k = std::floor(lz);

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};

    for (int ii = 0; ii <= 1; ++ii) { 
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int kk = 0; kk <= 1; ++kk) {
                amrex::Gpu::Atomic::Add(&rho(i+ii, j+jj, k+kk, 0),
                                        sx[ii]*sy[jj]*sz[kk]*p.m_rdata[0];
                for (int comp=1; comp < nc; ++comp) {
                    amrex::Gpu::Atomic::Add(&rho(i+ii, j+jj, k+kk, comp),
                                            sx[ii]*sy[jj]*sz[kk]*p.m_rdata[0]*p.m_rdata[comp];
                }
            }
        }
    }   
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
template <typename P>
void amrex_deposit_particle_dxcic (P const& p, int nc, amrex::Array4<amrex::Real const> const& rho,
                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                   amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& pdxi)
{
    amrex::Real factor = (pdxi[0]/dxi[0])*(pdxi[1]/dxi[1])*(pdxi[2]/dxi[2]);
    
    
    amrex::Abort("Not Implemented.");
}

AMREX_GPU_HOST_DEVICE AMREX_INLINE
template <typename P>
void amrex_interpolate_cic (P const& p, int nc, amrex::Array4<amrex::Real const> const& acc,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;
    
    int i = std::floor(lx);
    int j = std::floor(ly);
    int k = std::floor(lz);
    
    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;
    
    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};
    amrex::Real sz[] = {1.-zint, zint};
    
    for (int ii = 0; ii <= 1; ++ii) { 
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int kk = 0; kk <= 1; ++kk) {
                for (int comp=0; comp < nc; ++comp) {
                    acceleration[comp] = sx[ii]*sy[jj]*sz[kk]*acc(i+ii,j+jj,k+kk);
                }
            }
        }
    }
}

#endif
