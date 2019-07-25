#ifndef WARPX_GPU_PARSER_H_
#define WARPX_GPU_PARSER_H_

#include <WarpXParser.H>
#include <AMReX_Gpu.H>

class GpuParser
{
public:
    GpuParser (WarpXParser const& wp);
    void clear ();

    AMREX_GPU_HOST_DEVICE
    double
    operator() (double x, double y, double z) const noexcept
    {
#ifdef AMREX_USE_GPU

#ifdef AMREX_DEVICE_COMPILE
        amrex::Gpu::SharedMemory<double> gsm;
        double* p = gsm.dataPtr();
        int tid = threadIdx.x + threadIdx.y*blockDim.x + threadIdx.z*(blockDim.x*blockDim.y);
        p[tid*3] = x;
        p[tid*3+1] = y;
        p[tid*3+2] = z;
        return wp_ast_eval(m_gpu_parser.ast);
#else
        m_var.x = x;
        m_var.y = y;
        m_var.z = z;
        return wp_ast_eval(m_cpu_parser.ast);
#endif

#else

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
    struct wp_parser m_gpu_parser;
    struct wp_parser m_cpu_parser;
    mutable amrex::XDim3 m_var;
#else
    struct wp_parser** m_parser;
    mutable amrex::XDim3* m_var;
    int nthreads;
#endif
};

#endif
