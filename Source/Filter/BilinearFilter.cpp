#include <WarpX.H>
#include <BilinearFilter.H>
#include <WarpX_f.H>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace amrex;

namespace {
    void compute_stencil(Gpu::ManagedVector<Real> &stencil, int npass)
    {
        Gpu::ManagedVector<Real> old_s(1+npass,0.);
        Gpu::ManagedVector<Real> new_s(1+npass,0.);

        old_s[0] = 1.;
        int jmax = 1;
        amrex::Real loc;
        // Convolve the filter with itself npass times
        for(int ipass=1; ipass<npass+1; ipass++){
            // element 0 has to be treated in its own way
            new_s[0] = 0.5 * old_s[0];
            if (1<jmax) new_s[0] += 0.5 * old_s[1];
            loc = 0.;
            // For each element j, apply the filter to 
            // old_s to get new_s[j]. loc stores the tmp 
            // filtered value.
            for(int j=1; j<jmax+1; j++){
                loc = 0.5 * old_s[j];
                loc += 0.25 * old_s[j-1];
                if (j<jmax) loc += 0.25 * old_s[j+1];
                new_s[j] = loc;
            }
            // copy new_s into old_s
            old_s = new_s;
            // extend the stencil length for next iteration
            jmax += 1;
        }
        // we use old_s here to make sure the stencil
        // is corrent even when npass = 0
        stencil = old_s;
        stencil[0] *= 0.5; // because we will use it twice
    }
}

void BilinearFilter::ComputeStencils(){
    BL_PROFILE("BilinearFilter::ComputeStencils()");
    stencil_length_each_dir = npass_each_dir;
    stencil_length_each_dir += 1.;
#if (AMREX_SPACEDIM == 3)
    // npass_each_dir = npass_x npass_y npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_y.resize( 1 + npass_each_dir[1] );
    stencil_z.resize( 1 + npass_each_dir[2] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_y, npass_each_dir[1]);
    compute_stencil(stencil_z, npass_each_dir[2]);
#elif (AMREX_SPACEDIM == 2)
    // npass_each_dir = npass_x npass_z
    stencil_x.resize( 1 + npass_each_dir[0] );
    stencil_z.resize( 1 + npass_each_dir[1] );
    compute_stencil(stencil_x, npass_each_dir[0]);
    compute_stencil(stencil_z, npass_each_dir[1]);
#endif
    slen = stencil_length_each_dir.dim3();
#if (AMREX_SPACEDIM == 2)
    slen.z = 1;
#endif
}
