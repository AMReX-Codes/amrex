#ifndef WARPX_GPU_PARSER_H_
#define WARPX_GPU_PARSER_H_

#include <WarpXParser.H>
#include <AMReX_Gpu.H>

// When compiled for CPU, wrap WarpXParser and enable threading.
// When compiled for GPU, store one copy of the parser in
// CUDA managed memory for __device__ code, and one copy of the parser
// in CUDA managed memory for __host__ code. This way, the parser can be
// efficiently called from both host and device.
class GpuParser
{
public:
    GpuParser (WarpXParser const& wp);
    void clear ();

    AMREX_GPU_HOST_DEVICE
    amrex::Real
    operator() (amrex::Real x, amrex::Real y, amrex::Real z) const noexcept
    {
#ifdef AMREX_USE_GPU

#ifdef AMREX_DEVICE_COMPILE
// WarpX compiled for GPU, function compiled for __device__
        // the 3D position of each particle is stored in shared memory.
        amrex::Gpu::SharedMemory<amrex::Real> gsm;
        amrex::Real* p = gsm.dataPtr();
        int tid = threadIdx.x + threadIdx.y*blockDim.x + threadIdx.z*(blockDim.x*blockDim.y);
        p[tid*3] = x;
        p[tid*3+1] = y;
        p[tid*3+2] = z;
        return wp_ast_eval(m_gpu_parser.ast);
#else
// WarpX compiled for GPU, function compiled for __host__
        m_var.x = x;
        m_var.y = y;
        m_var.z = z;
        return wp_ast_eval(m_cpu_parser.ast);
#endif

#else
// WarpX compiled for CPU
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        m_var[tid].x = x;
        m_var[tid].y = y;
        m_var[tid].z = z;
        return wp_ast_eval(m_parser[tid]->ast);
#endif
    }

private:

#ifdef AMREX_USE_GPU
    // Copy of the parser running on __device__
    struct wp_parser m_gpu_parser;
    // Copy of the parser running on __host__
    struct wp_parser m_cpu_parser;
    mutable amrex::XDim3 m_var;
#else
    // Only one parser
    struct wp_parser** m_parser;
    mutable amrex::XDim3* m_var;
    int nthreads;
#endif
};

#endif
