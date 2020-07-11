#include <AMReX_Gpu.H>

namespace amrex {
namespace Gpu {

#ifdef AMREX_USE_CUDA

Fuser3D::Fuser3D ()
{
    m_lambda_buf = (char*)The_Pinned_Arena()->alloc(m_nbytes_lambda_buf);
    m_helper_buf = (Fuse3DHelper*)The_Pinned_Arena()->alloc(m_nhelpers_buf*sizeof(Fuse3DHelper));
}

Fuser3D::~Fuser3D ()
{
    if (m_nlambdas > 0){
        Launch();
    }

    The_Pinned_Arena()->free(m_lambda_buf);
    The_Pinned_Arena()->free(m_helper_buf);
}

void Fuser3D::Launch ()
{
    int nlambdas = m_nlambdas;
    if (nlambdas > 0) {
        int* nwarps = (int*)The_Pinned_Arena()->alloc(nlambdas*sizeof(int));
        int ntotwarps = 0;
        for (int i = 0; i < nlambdas; ++i)
        {
            nwarps[i] = ntotwarps;
            ntotwarps += static_cast<int>(m_helper_buf[i].m_bx.numPts()+Gpu::Device::warp_size-1)
                / Gpu::Device::warp_size;
        }
        nwarps[nlambdas] = ntotwarps;

        int* d_nwarps = (int*)The_Device_Arena()->alloc((nlambdas+1)*sizeof(int));
        auto d_lambda_helper = (Fuse3DHelper*)The_Device_Arena()->alloc
            (nlambdas*sizeof(Fuse3DHelper));
        auto d_lambda_object = (char*)The_Device_Arena()->alloc(m_nbytes_used_lambda_buf);

        Gpu::htod_memcpy_async(d_nwarps, nwarps, (nlambdas+1)*sizeof(int));
        Gpu::htod_memcpy_async(d_lambda_helper, m_helper_buf, nlambdas*sizeof(Fuse3DHelper));
        Gpu::htod_memcpy_async(d_lambda_object, m_lambda_buf, m_nbytes_used_lambda_buf);

        constexpr int nthreads = 256;
        constexpr int nwarps_per_block = nthreads/Gpu::Device::warp_size;
        int nblocks = (ntotwarps + nwarps_per_block-1) / nwarps_per_block;

        amrex::launch(nblocks, nthreads, Gpu::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int g_tid = blockDim.x*blockIdx.x + threadIdx.x;
            int g_wid = g_tid / Gpu::Device::warp_size;
            if (g_wid >= ntotwarps) return;

            int ilambda;
            {
                int lo = 0;
                int hi = nlambdas;
                while (lo <= hi) {
                    int mid = (lo+hi)/2;
                    if (g_wid >= d_nwarps[mid] && g_wid < d_nwarps[mid+1]) {
                        ilambda = mid;
                        break;
                    } else if (g_wid < d_nwarps[mid]) {
                        hi = mid-1;
                    } else {
                        lo = mid+1;
                    }
                };
            }

            Box const& bx = d_lambda_helper[ilambda].m_bx;
            int ncells = bx.numPts();
            int b_wid = g_wid - d_nwarps[ilambda]; // b_wid'th warp on this this lambda
            int lane = threadIdx.x % Gpu::Device::warp_size;
            int icell = b_wid*Gpu::Device::warp_size + lane;
            if (icell < ncells) {
                const auto len = amrex::length(bx);
                const auto lo  = amrex::lbound(bx);
                int k =  icell /   (len.x*len.y);
                int j = (icell - k*(len.x*len.y)) /   len.x;
                int i = (icell - k*(len.x*len.y)) - j*len.x;
                i += lo.x;
                j += lo.y;
                k += lo.z;
                d_lambda_helper[ilambda].m_fp(d_lambda_object+d_lambda_helper[ilambda].m_offset,
                                              i,j,k);
            }
        });
        Gpu::synchronize();
        The_Pinned_Arena()->free(nwarps);
        The_Device_Arena()->free(d_nwarps);
        The_Device_Arena()->free(d_lambda_helper);
        The_Device_Arena()->free(d_lambda_object);
        m_nlambdas = 0;
    }

}

#endif

}}
