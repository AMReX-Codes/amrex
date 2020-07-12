#include <AMReX_Gpu.H>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {
namespace Gpu {

#ifdef AMREX_USE_CUDA

Fuser::Fuser ()
{
    m_lambda_buf = (char*)The_Pinned_Arena()->alloc(m_nbytes_lambda_buf);
    m_helper_buf = (FuseHelper*)The_Pinned_Arena()->alloc(m_nhelpers_buf*sizeof(FuseHelper));
    m_dtor_buf.reserve(1024);
#ifdef _OPENMP
    AMREX_ASSERT(!omp_in_parallel());
#endif
}

Fuser::~Fuser ()
{
    if (m_nlambdas > 0){
        Launch();
    }

    The_Pinned_Arena()->free(m_lambda_buf);
    The_Pinned_Arena()->free(m_helper_buf);
}

void Fuser::Launch ()
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
        auto d_lambda_helper = (FuseHelper*)The_Device_Arena()->alloc
            (nlambdas*sizeof(FuseHelper));
        auto d_lambda_object = (char*)The_Device_Arena()->alloc(m_nbytes_used_lambda_buf);

        Gpu::htod_memcpy_async(d_nwarps, nwarps, (nlambdas+1)*sizeof(int));
        Gpu::htod_memcpy_async(d_lambda_helper, m_helper_buf, nlambdas*sizeof(FuseHelper));
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
                FuseHelper& helper = d_lambda_helper[ilambda];
                if (helper.m_ncomp == 0) {
                    helper.m_fp.L3D(d_lambda_object+helper.m_offset,i,j,k);
                } else {
                    for (int n = 0; n < helper.m_ncomp; ++n) {
                        helper.m_fp.L4D(d_lambda_object+helper.m_offset,i,j,k,n);
                    }
                }
            }
        });
        Gpu::synchronize();
        The_Pinned_Arena()->free(nwarps);
        The_Device_Arena()->free(d_nwarps);
        The_Device_Arena()->free(d_lambda_helper);
        The_Device_Arena()->free(d_lambda_object);

        for (int i = 0; i < nlambdas; ++i) {
            char* p = m_lambda_buf + m_helper_buf[i].m_offset;
            m_dtor_buf[i](p);
            m_helper_buf[i].~FuseHelper();
        }
        m_dtor_buf.clear();
        m_nlambdas = 0;
    }
}

void
Fuser::resize_lambda_buf ()
{
    m_nbytes_lambda_buf += m_nbytes_lambda_buf/2;
    auto p = (char*)The_Pinned_Arena()->alloc(m_nbytes_lambda_buf);
    std::memcpy(p, m_lambda_buf, m_nbytes_used_lambda_buf);
    The_Pinned_Arena()->free(m_lambda_buf);
    m_lambda_buf = p;
}

void
Fuser::resize_helper_buf ()
{
    m_nhelpers_buf += m_nhelpers_buf/2;
    auto p = (FuseHelper*)The_Pinned_Arena()->alloc(m_nhelpers_buf*sizeof(FuseHelper));
    for (int i = 0; i < m_nlambdas; ++i) {
        new (p+i) FuseHelper(m_helper_buf[i]);
        (m_helper_buf+i)->~FuseHelper();
    }
    The_Pinned_Arena()->free(m_helper_buf);
    m_helper_buf = p;
}

#endif

}}
