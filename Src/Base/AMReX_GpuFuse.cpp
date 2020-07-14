#include <AMReX_Gpu.H>
#include <AMReX_ParmParse.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_OpenMP.H>

namespace amrex {
namespace Gpu {

#ifdef AMREX_USE_CUDA

namespace {
    bool s_in_fuse_region = false;
    // Fusing kernels with elements greater than this are not recommended based tests on v100
    Long s_fuse_size_threshold = 257*257;
    // If the number of kernels is less than this, fusing is not recommended based on tests on v100
    int s_fuse_numkernels_threshold = 3;
}

std::unique_ptr<Fuser> Fuser::m_instance = nullptr;

Fuser::Fuser ()
{
    AMREX_ASSERT(!OpenMP::in_parallel());
    m_lambda_buf = (char*)The_Pinned_Arena()->alloc(m_nbytes_lambda_buf);
    m_helper_buf = (FuseHelper*)The_Pinned_Arena()->alloc(m_nhelpers_buf*sizeof(FuseHelper));
    m_dtor_buf.reserve(1024);
}

Fuser::~Fuser ()
{
    if (m_nlambdas > 0) {
        Launch();
    }

    The_Pinned_Arena()->free(m_lambda_buf);
    The_Pinned_Arena()->free(m_helper_buf);
}

void Fuser::Launch ()
{
    BL_PROFILE("Fuser::Launch()");

    int nlambdas = m_nlambdas;
    if (nlambdas > 0) {
        AMREX_ASSERT(!OpenMP::in_parallel());

        int* nwarps = (int*)The_Pinned_Arena()->alloc((nlambdas+1)*sizeof(int));
        int ntotwarps = 0;
        for (int i = 0; i < nlambdas; ++i)
        {
            nwarps[i] = ntotwarps;
            Box const& bx = m_helper_buf[i].m_bx;
            int N;
            if (bx.isEmpty()) {
                N = m_helper_buf[i].m_N;
            } else {
                N = bx.numPts();
            }
            ntotwarps += (N + Gpu::Device::warp_size-1) / Gpu::Device::warp_size;
        }
        nwarps[nlambdas] = ntotwarps;

        AMREX_ASSERT(ntotwarps < std::numeric_limits<int>::max()/Gpu::Device::warp_size);

        // Pack nwarps, lambda helpers and lambda objects into a buffer
        std::size_t sizeof_nwarps = sizeof(int) * (nlambdas+1);
        std::size_t offset_helpers = Arena::align(sizeof_nwarps);
        std::size_t sizeof_helpers = sizeof(FuseHelper)*nlambdas;
        std::size_t offset_objects = Arena::align(offset_helpers+sizeof_helpers);
        std::size_t sizeof_objects = m_nbytes_used_lambda_buf;
        std::size_t total_buf_size = offset_objects + sizeof_objects;

        char* h_buffer = (char*)The_Pinned_Arena()->alloc(total_buf_size);
        char* d_buffer = (char*)The_Device_Arena()->alloc(total_buf_size);

        std::memcpy(h_buffer, nwarps, sizeof_nwarps);
        std::memcpy(h_buffer+offset_helpers, m_helper_buf, sizeof_helpers);
        std::memcpy(h_buffer+offset_objects, m_lambda_buf, sizeof_objects);
        Gpu::htod_memcpy_async(d_buffer, h_buffer, total_buf_size);

        auto d_nwarps = reinterpret_cast<int*>(d_buffer);
        auto d_lambda_helper = reinterpret_cast<FuseHelper*>(d_buffer+offset_helpers);
        auto d_lambda_object = reinterpret_cast<char*>(d_buffer+offset_objects);

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

            int b_wid = g_wid - d_nwarps[ilambda]; // b_wid'th warp on this lambda
            int lane = threadIdx.x % Gpu::Device::warp_size;
            int icell = b_wid*Gpu::Device::warp_size + lane;

            FuseHelper& helper = d_lambda_helper[ilambda];
            char* lambda_object = d_lambda_object + helper.m_offset;
            Box const& bx = helper.m_bx;
            if (bx.isEmpty()) {
                if (icell < helper.m_N) {
                    helper.m_fp.L1D(lambda_object,icell);
                }
            } else {
                int ncells = bx.numPts();
                if (icell < ncells) {
                    const auto len = amrex::length(bx);
                    const auto lo  = amrex::lbound(bx);
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    if (helper.m_N == 0) {
                        helper.m_fp.L3D(lambda_object,i,j,k);
                    } else {
                        for (int n = 0; n < helper.m_N; ++n) {
                            helper.m_fp.L4D(lambda_object,i,j,k,n);
                        }
                    }
                }
            }
        });
        Gpu::synchronize();
        The_Pinned_Arena()->free(nwarps);
        The_Pinned_Arena()->free(h_buffer);
        The_Device_Arena()->free(d_buffer);

        for (int i = 0; i < nlambdas; ++i) {
            char* p = m_lambda_buf + m_helper_buf[i].m_offset;
            m_dtor_buf[i](p);
            m_helper_buf[i].~FuseHelper();
        }
        m_dtor_buf.clear();
        m_nbytes_used_lambda_buf = 0;
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

Fuser&
Fuser::getInstance ()
{
    if (m_instance == nullptr) {
        m_instance.reset(new Fuser());
    }
    return *m_instance;
}

void
Fuser::Initialize ()
{
    ParmParse pp("amrex");
    pp.query("gpu_fuse_size_threshold", s_fuse_size_threshold);
    pp.query("gpu_fuse_numkernels_threshold", s_fuse_numkernels_threshold);

    amrex::ExecOnFinalize(Fuser::Finalize);
}

void
Fuser::Finalize ()
{
    m_instance.reset();
}

Long getFuseSizeThreshold () { return s_fuse_size_threshold; }

Long
setFuseSizeThreshold (Long new_threshold)
{
    Long old = s_fuse_size_threshold;
    s_fuse_size_threshold = new_threshold;
    return old;
}

int getFuseNumKernelsThreshold () { return s_fuse_numkernels_threshold; }

int
setFuseNumKernelsThreshold (int new_threshold)
{
    int old = s_fuse_numkernels_threshold;
    s_fuse_numkernels_threshold = new_threshold;
    return old;
}

bool inFuseRegion () { return s_in_fuse_region; }

bool
setFuseRegion (bool flag)
{
    bool old = s_in_fuse_region;
    s_in_fuse_region = flag;
    return old;
}

#endif

}}
