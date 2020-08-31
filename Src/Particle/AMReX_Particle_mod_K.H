#ifndef AMREX_PARTICLE_MOD_K_H_
#define AMREX_PARTICLE_MOD_K_H_

#include <AMReX_Particles.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_REAL.H>
#include <cmath>

namespace amrex {

template <typename P>
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void amrex_deposit_cic (P const& p, int nc, amrex::Array4<amrex::Real> const& rho,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& plo,
                        amrex::GpuArray<amrex::Real,AMREX_SPACEDIM> const& dxi)
{
#if (AMREX_SPACEDIM == 1)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    
    int i = static_cast<int>(amrex::Math::floor(lx));
    
    amrex::Real xint = lx - i;
    
    amrex::Real sx[2] = {Real(1.0) - xint, xint};
    
    for (int ii = 0; ii <= 1; ++ii) { 
        amrex::Gpu::Atomic::Add(&rho(i+ii-1, 0, 0, 0), static_cast<Real>(sx[ii]*p.rdata(0)));
    }
    
    for (int comp=1; comp < nc; ++comp) {
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Gpu::Atomic::Add(&rho(i+ii-1, 0, 0, comp), 
                                    static_cast<Real>(sx[ii]*p.rdata(0)*p.rdata(comp)));
        }
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;

    amrex::Real sx[2] = {Real(1.0) - xint, xint};
    amrex::Real sy[2] = {Real(1.0) - yint, yint};

    for (int jj = 0; jj <= 1; ++jj) { 
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, 0, 0),
                                    static_cast<Real>(sx[ii]*sy[jj]*p.rdata(0)));
        }
    }

    for (int comp=1; comp < nc; ++comp) {
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) {                 
                amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, 0, comp),
                                        static_cast<Real>(sx[ii]*sy[jj]*p.rdata(0)*p.rdata(comp)));
            }
        }
    }

#elif (AMREX_SPACEDIM == 3)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;

    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));

    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;

    amrex::Real sx[] = {Real(1.0) - xint, xint};
    amrex::Real sy[] = {Real(1.0) - yint, yint};
    amrex::Real sz[] = {Real(1.0) - zint, zint};
    
    for (int kk = 0; kk <= 1; ++kk) { 
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) {
                amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, 0),
                                        static_cast<Real>(sx[ii]*sy[jj]*sz[kk]*p.rdata(0)));
            }
        }
    }

    for (int comp=1; comp < nc; ++comp) {
        for (int kk = 0; kk <= 1; ++kk) { 
            for (int jj = 0; jj <= 1; ++jj) { 
                for (int ii = 0; ii <= 1; ++ii) {
                    amrex::Gpu::Atomic::Add(&rho(i+ii-1, j+jj-1, k+kk-1, comp),
                                            static_cast<Real>(sx[ii]*sy[jj]*sz[kk]*p.rdata(0)*p.rdata(comp)));
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

    int lo_x = static_cast<int>(amrex::Math::floor(lx));

    int hi_x = static_cast<int>(amrex::Math::floor(hx));
    
    for (int i = lo_x; i <= hi_x; ++i) {
        if (i < rho.begin.x || i >= rho.end.x) continue;
        amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
        amrex::Real weight = wx*factor;                
        rho(i, 0, 0, 0) += weight*p.rdata(0);                
    }

    for (int comp = 1; comp < nc; ++comp)
    {
        for (int i = lo_x; i <= hi_x; ++i) {
            if (i < rho.begin.x || i >= rho.end.x) continue;
            amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));            
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

    int lo_x = static_cast<int>(amrex::Math::floor(lx));
    int lo_y = static_cast<int>(amrex::Math::floor(ly));

    int hi_x = static_cast<int>(amrex::Math::floor(hx));
    int hi_y = static_cast<int>(amrex::Math::floor(hy));

    for (int j = lo_y; j <= hi_y; ++j) {
        if (j < rho.begin.y || j >= rho.end.y) continue;
        amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
        for (int i = lo_x; i <= hi_x; ++i) {
            if (i < rho.begin.x || i >= rho.end.x) continue;
            amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
            
            amrex::Real weight = wx*wy*factor;
                
            rho(i, j, 0, 0) += weight*p.rdata(0);                
        }
    }

    for (int comp = 1; comp < nc; ++comp) {
        for (int j = lo_y; j <= hi_y; ++j) {
            if (j < rho.begin.y || j >= rho.end.y) continue;
            amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
            for (int i = lo_x; i <= hi_x; ++i) {
                if (i < rho.begin.x || i >= rho.end.x) continue;
                amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
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
    
    int lo_x = static_cast<int>(amrex::Math::floor(lx));
    int lo_y = static_cast<int>(amrex::Math::floor(ly));
    int lo_z = static_cast<int>(amrex::Math::floor(lz));
    
    int hi_x = static_cast<int>(amrex::Math::floor(hx));
    int hi_y = static_cast<int>(amrex::Math::floor(hy));
    int hi_z = static_cast<int>(amrex::Math::floor(hz));
    
    for (int k = lo_z; k <= hi_z; ++k) {
        if (k < rho.begin.z || k >= rho.end.z) continue;
        amrex::Real wz = amrex::min(hz - k, amrex::Real(1.0)) - amrex::max(lz - k, amrex::Real(0.0));
        for (int j = lo_y; j <= hi_y; ++j) {
            if (j < rho.begin.y || j >= rho.end.y) continue;
            amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
            for (int i = lo_x; i <= hi_x; ++i) {
                if (i < rho.begin.x || i >= rho.end.x) continue;
                amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
                
                amrex::Real weight = wx*wy*wz*factor;
                
                rho(i, j, k, 0) += weight*p.rdata(0);                
            }
        }
    }

    for (int comp = 1; comp < nc; ++comp) {
        for (int k = lo_z; k <= hi_z; ++k) {
            if (k < rho.begin.z || k >= rho.end.z) continue;
            amrex::Real wz = amrex::min(hz - k, amrex::Real(1.0)) - amrex::max(lz - k, amrex::Real(0.0));
            for (int j = lo_y; j <= hi_y; ++j) {
                if (j < rho.begin.y || j >= rho.end.y) continue;
                amrex::Real wy = amrex::min(hy - j, amrex::Real(1.0)) - amrex::max(ly - j, amrex::Real(0.0));
                for (int i = lo_x; i <= hi_x; ++i) {
                    if (i < rho.begin.x || i >= rho.end.x) continue;
                    amrex::Real wx = amrex::min(hx - i, amrex::Real(1.0)) - amrex::max(lx - i, amrex::Real(0.0));
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
    
    int i = static_cast<int>(amrex::Math::floor(lx));
    
    amrex::Real xint = lx - i;
    
    amrex::Real sx[] = {Real(1.0)-xint, xint};

    for (int comp=0; comp < nc; ++comp) {    
        for (int ii = 0; ii <= 1; ++ii) { 
            amrex::Real acceleration = sx[ii]*acc(i+ii-1,0,0,comp);
            amrex::ignore_unused(acceleration);
        }
    }
#elif (AMREX_SPACEDIM == 2)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    
    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    
    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    
    amrex::Real sx[] = {Real(1.0)-xint, xint};
    amrex::Real sy[] = {Real(1.0)-yint, yint};

    for (int comp=0; comp < nc; ++comp) {    
        for (int jj = 0; jj <= 1; ++jj) { 
            for (int ii = 0; ii <= 1; ++ii) { 
                amrex::Real acceleration = sx[ii]*sy[jj]*acc(i+ii-1,j+jj-1,0,comp);
                amrex::ignore_unused(acceleration);
            }
        }
    }
#elif (AMREX_SPACEDIM == 3)
    amrex::Real lx = (p.pos(0) - plo[0]) * dxi[0] + 0.5;
    amrex::Real ly = (p.pos(1) - plo[1]) * dxi[1] + 0.5;
    amrex::Real lz = (p.pos(2) - plo[2]) * dxi[2] + 0.5;
    
    int i = static_cast<int>(amrex::Math::floor(lx));
    int j = static_cast<int>(amrex::Math::floor(ly));
    int k = static_cast<int>(amrex::Math::floor(lz));
    
    amrex::Real xint = lx - i;
    amrex::Real yint = ly - j;
    amrex::Real zint = lz - k;
    
    amrex::Real sx[] = {Real(1.0)-xint, xint};
    amrex::Real sy[] = {Real(1.0)-yint, yint};
    amrex::Real sz[] = {Real(1.0)-zint, zint};

    for (int comp=0; comp < nc; ++comp) {                    
        for (int kk = 0; kk <= 1; ++kk) { 
            for (int jj = 0; jj <= 1; ++jj) { 
                for (int ii = 0; ii <= 1; ++ii) {
                    amrex::Real acceleration = sx[ii]*sy[jj]*sz[kk]*acc(i+ii-1,j+jj-1,k+kk-1,comp);
                    amrex::ignore_unused(acceleration);
                }
            }
        }
    }
#else
    amrex::Abort("Not implemented.");
#endif
}

}

#endif
