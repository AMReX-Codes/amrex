#ifndef AMREX_PARTICLE_MOD_K_H_
#define AMREX_PARTICLE_MOD_K_H_

#include <AMReX_Particles.H>
#include <AMReX_FArrayBox.H>
#include <cmath>

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_deposit_cic (P const& p, int nc, amrex::Array4<amrex::Real> const& rho,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
#if (AMREX_SPACEDIM == 1)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    
    int i = std::floor(lx);
    
    amrex::Real xint = lx - i;
    
    amrex::Real sx[] = {1.-xint, xint};
    
    for (int ii = 0; ii <= 1; ++ii) { 
        amrex::Gpu::Atomic::Add(&rho(i+ii-1, 0, 0, 0), sx[ii]*p.rdata(0));
    }
    
    for (int comp=1; comp < nc; ++comp) {
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Gpu::Atomic::Add(&rho(i+ii-1, 0, 0, comp), sx[ii]*p.rdata(0)*p.rdata(comp));
        }
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;

    int i = std::floor(lx);
    int j = std::floor(ly);

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;

    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};

    for (int jj = 0; jj <= 1; ++jj) { 
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, 0, 0),
                                    sx[ii]*sy[jj]*p.rdata(0));
        }
    }

    for (int comp=1; comp < nc; ++comp) {
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) {                 
                amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, 0, comp),
                                        sx[ii]*sy[jj]*p.rdata(0)*p.rdata(comp));
            }
        }
    }

#elif (AMREX_SPACEDIM == 3)
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
    
    for (int kk = 0; kk <= 1; ++kk) { 
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) {
                amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, 0),
                                        sx[ii]*sy[jj]*sz[kk]*p.rdata(0));
            }
        }
    }

    for (int comp=1; comp < nc; ++comp) {
        for (int kk = 0; kk <= 1; ++kk) { 
            for (int jj = 0; jj <= 1; ++jj) { 
                for (int ii = 0; ii <= 1; ++ii) {
                    amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, comp),
                                            sx[ii]*sy[jj]*sz[kk]*p.rdata(0)*p.rdata(comp));
                }
            }
        }
    }
#else
    amrex::Abort("Not implemented.");
#endif
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_deposit_particle_dx_cic (P const& p, int nc, amrex::Array4<amrex::Real> const& rho,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi,
                                    amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& pdxi)
{
#if (AMREX_SPACEDIM == 1)
    amrex::Real factor = (pdxi[0]/dxi[0]);

    amrex::Real lx = (p.pos(0) - plo[0] - 0.5/pdxi[0]) * dxi[0];

    amrex::Real hx = (p.pos(0) - plo[0] + 0.5/pdxi[0]) * dxi[0];

    int lo_x = std::floor(lx);

    int hi_x = std::floor(hx);
    
    for (int i = lo_x; i <= hi_x; ++i) {
        if (i < rho.begin.x || i >= rho.end.x) continue;
        amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);                
        amrex::Real weight = wx*factor;                
        rho(i, 0, 0, 0) += weight*p.rdata(0);                
    }

    for (int comp = 1; comp < nc; ++comp)
    {
        for (int i = lo_x; i <= hi_x; ++i) {
            if (i < rho.begin.x || i >= rho.end.x) continue;
            amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);            
            amrex::Real weight = wx*factor;
            rho(i, 0, 0, comp) += weight*p.rdata(0)*p.rdata(comp);
        }
    }

#elif (AMREX_SPACEDIM == 2)
    amrex::Real factor = (pdxi[0]/dxi[0])*(pdxi[1]/dxi[1]);

    amrex::Real lx = (p.pos(0) - plo[0] - 0.5/pdxi[0]) * dxi[0];
    amrex::Real ly = (p.pos(1) - plo[1] - 0.5/pdxi[1]) * dxi[1];

    amrex::Real hx = (p.pos(0) - plo[0] + 0.5/pdxi[0]) * dxi[0];
    amrex::Real hy = (p.pos(1) - plo[1] + 0.5/pdxi[1]) * dxi[1];

    int lo_x = std::floor(lx);
    int lo_y = std::floor(ly);

    int hi_x = std::floor(hx);
    int hi_y = std::floor(hy);

    for (int j = lo_y; j <= hi_y; ++j) {
        if (j < rho.begin.y || j >= rho.end.y) continue;
        amrex::Real wy = amrex::min(hy - j, 1.0) - amrex::max(ly - j, 0.0);
        for (int i = lo_x; i <= hi_x; ++i) {
            if (i < rho.begin.x || i >= rho.end.x) continue;
            amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);
            
            amrex::Real weight = wx*wy*factor;
                
            rho(i, j, 0, 0) += weight*p.rdata(0);                
        }
    }

    for (int comp = 1; comp < nc; ++comp) {
        for (int j = lo_y; j <= hi_y; ++j) {
            if (j < rho.begin.y || j >= rho.end.y) continue;
            amrex::Real wy = amrex::min(hy - j, 1.0) - amrex::max(ly - j, 0.0);
            for (int i = lo_x; i <= hi_x; ++i) {
                if (i < rho.begin.x || i >= rho.end.x) continue;
                amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);
                amrex::Real weight = wx*wy*factor;
                rho(i, j, 0, comp) += weight*p.rdata(0)*p.rdata(comp);
            }
        }
    }

#elif (AMREX_SPACEDIM == 3)
    amrex::Real factor = (pdxi[0]/dxi[0])*(pdxi[1]/dxi[1])*(pdxi[2]/dxi[2]);
    
    amrex::Real lx = (p.pos(0) - plo[0] - 0.5/pdxi[0]) * dxi[0];
    amrex::Real ly = (p.pos(1) - plo[1] - 0.5/pdxi[1]) * dxi[1];
    amrex::Real lz = (p.pos(2) - plo[2] - 0.5/pdxi[2]) * dxi[2];
    
    amrex::Real hx = (p.pos(0) - plo[0] + 0.5/pdxi[0]) * dxi[0];
    amrex::Real hy = (p.pos(1) - plo[1] + 0.5/pdxi[1]) * dxi[1];
    amrex::Real hz = (p.pos(2) - plo[2] + 0.5/pdxi[2]) * dxi[2];
    
    int lo_x = std::floor(lx);
    int lo_y = std::floor(ly);
    int lo_z = std::floor(lz);
    
    int hi_x = std::floor(hx);
    int hi_y = std::floor(hy);
    int hi_z = std::floor(hz);
    
    for (int k = lo_z; k <= hi_z; ++k) {
        if (k < rho.begin.z || k >= rho.end.z) continue;
        amrex::Real wz = amrex::min(hz - k, 1.0) - amrex::max(lz - k, 0.0);
        for (int j = lo_y; j <= hi_y; ++j) {
            if (j < rho.begin.y || j >= rho.end.y) continue;
            amrex::Real wy = amrex::min(hy - j, 1.0) - amrex::max(ly - j, 0.0);
            for (int i = lo_x; i <= hi_x; ++i) {
                if (i < rho.begin.x || i >= rho.end.x) continue;
                amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);
                
                amrex::Real weight = wx*wy*wz*factor;
                
                rho(i, j, k, 0) += weight*p.rdata(0);                
            }
        }
    }

    for (int comp = 1; comp < nc; ++comp) {
        for (int k = lo_z; k <= hi_z; ++k) {
            if (k < rho.begin.z || k >= rho.end.z) continue;
            amrex::Real wz = amrex::min(hz - k, 1.0) - amrex::max(lz - k, 0.0);
            for (int j = lo_y; j <= hi_y; ++j) {
                if (j < rho.begin.y || j >= rho.end.y) continue;
                amrex::Real wy = amrex::min(hy - j, 1.0) - amrex::max(ly - j, 0.0);
                for (int i = lo_x; i <= hi_x; ++i) {
                    if (i < rho.begin.x || i >= rho.end.x) continue;
                    amrex::Real wx = amrex::min(hx - i, 1.0) - amrex::max(lx - i, 0.0);
                    amrex::Real weight = wx*wy*wz*factor;
                    rho(i, j, k, comp) += weight*p.rdata(0)*p.rdata(comp);
                }
            }
        }
    }
#else
    amrex::Abort("Not implemented.")
#endif
}

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_interpolate_cic (P const& p, int nc, amrex::Array4<amrex::Real const> const& acc,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                            amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
#if (AMREX_SPACEDIM == 1)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    
    int i = std::floor(lx);
    
    amrex::Real xint = lx - i;
    
    amrex::Real sx[] = {1.-xint, xint};

    for (int comp=0; comp < nc; ++comp) {    
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Real acceleration = sx[ii]*acc(i+ii,0,0,comp);
        }
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    
    int i = std::floor(lx);
    int j = std::floor(ly);
    
    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    
    amrex::Real sx[] = {1.-xint, xint};
    amrex::Real sy[] = {1.-yint, yint};

    for (int comp=0; comp < nc; ++comp) {    
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) { 
                amrex::Real acceleration = sx[ii]*sy[jj]*acc(i+ii,j+jj,0,comp);
            }
        }
    }
#elif (AMREX_SPACEDIM == 3)
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

    for (int comp=0; comp < nc; ++comp) {                    
        for (int kk = 0; kk <= 1; ++kk) { 
            for (int jj = 0; jj <= 1; ++jj) { 
                for (int ii = 0; ii <= 1; ++ii) {
                    amrex::Real acceleration = sx[ii]*sy[jj]*sz[kk]*acc(i+ii,j+jj,k+kk,comp);
                }
            }
        }
    }
#else
    amrex::Abort("Not implemented.");
#endif
}

#endif
