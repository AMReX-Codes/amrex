/* Copyright 2019-2020 Maxence Thevenet, Revathi Jambunathan, Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_GPU_PARSER_H_
#define WARPX_GPU_PARSER_H_

#include "Parser/WarpXParser.H"

#include <AMReX_Gpu.H>
#include <AMReX_Array.H>
#include <AMReX_TypeTraits.H>


// When compiled for CPU, wrap WarpXParser and enable threading.
// When compiled for GPU, store one copy of the parser in
// CUDA managed memory for __device__ code, and one copy of the parser
// in CUDA managed memory for __host__ code. This way, the parser can be
// efficiently called from both host and device.
template <int N>
class GpuParser
{
public:
    GpuParser (WarpXParser const& wp);
    void clear ();

    template <typename... Ts>
    AMREX_GPU_HOST_DEVICE
    std::enable_if_t<sizeof...(Ts) == N
                     and amrex::Same<amrex::Real,Ts...>::value,
                     amrex::Real>
    operator() (Ts... var) const noexcept
    {
#ifdef AMREX_USE_GPU
        amrex::GpuArray<amrex::Real,N> l_var{var...};
#if defined(__CUDA_ARCH__) || defined(__HIP_DEVICE_COMPILE__)
// WarpX compiled for GPU, function compiled for __device__
        return wp_ast_eval(m_gpu_parser.ast, l_var.data());
#else
// WarpX compiled for GPU, function compiled for __host__
        return wp_ast_eval(m_cpu_parser->ast, nullptr);
#endif

#else
// WarpX compiled for CPU
#ifdef _OPENMP
        int tid = omp_get_thread_num();
#else
        int tid = 0;
#endif
        m_var[tid] = amrex::GpuArray<amrex::Real,N>{var...};
        return wp_ast_eval(m_parser[tid]->ast, nullptr);
#endif
    }


private:

#ifdef AMREX_USE_GPU
    // Copy of the parser running on __device__
    struct wp_parser m_gpu_parser;
    // Copy of the parser running on __host__
    struct wp_parser* m_cpu_parser;
    mutable amrex::GpuArray<amrex::Real,N> m_var;
#else
    // Only one parser
    struct wp_parser** m_parser;
    mutable amrex::GpuArray<amrex::Real,N>* m_var;
    int nthreads;
#endif
};

template <int N>
GpuParser<N>::GpuParser (WarpXParser const& wp)
{
#ifdef AMREX_USE_GPU

    struct wp_parser* a_wp = wp.m_parser;
    // Initialize GPU parser: allocate memory in CUDA managed memory,
    // copy all data needed on GPU to m_gpu_parser
    m_gpu_parser.sz_mempool = wp_ast_size(a_wp->ast);
    m_gpu_parser.p_root = (struct wp_node*)
        amrex::The_Managed_Arena()->alloc(m_gpu_parser.sz_mempool);
    m_gpu_parser.p_free = m_gpu_parser.p_root;
    // 0: don't free the source
    m_gpu_parser.ast = wp_parser_ast_dup(&m_gpu_parser, a_wp->ast, 0);
    for (int i = 0; i < N; ++i) {
        wp_parser_regvar_gpu(&m_gpu_parser, wp.m_varnames[i].c_str(), i);
    }

    // Initialize CPU parser:
    m_cpu_parser = wp_parser_dup(a_wp);
    for (int i = 0; i < N; ++i) {
        wp_parser_regvar(m_cpu_parser, wp.m_varnames[i].c_str(), &m_var[i]);
    }

#else // not defined AMREX_USE_GPU

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else // _OPENMP
    nthreads = 1;
#endif // _OPENMP

    m_parser = ::new struct wp_parser*[nthreads];
    m_var = ::new amrex::GpuArray<amrex::Real,N>[nthreads];

    for (int tid = 0; tid < nthreads; ++tid)
    {
#ifdef _OPENMP
        m_parser[tid] = wp_parser_dup(wp.m_parser[tid]);
        for (int i = 0; i < N; ++i) {
            wp_parser_regvar(m_parser[tid], wp.m_varnames[tid][i].c_str(), &(m_var[tid][i]));
        }
#else // _OPENMP
        m_parser[tid] = wp_parser_dup(wp.m_parser);
        for (int i = 0; i < N; ++i) {
            wp_parser_regvar(m_parser[tid], wp.m_varnames[i].c_str(), &(m_var[tid][i]));
        }
#endif // _OPENMP
    }

#endif // AMREX_USE_GPU
}


template <int N>
void
GpuParser<N>::clear ()
{
#ifdef AMREX_USE_GPU
    amrex::The_Managed_Arena()->free(m_gpu_parser.ast);
    wp_parser_delete(m_cpu_parser);
#else
    for (int tid = 0; tid < nthreads; ++tid)
    {
        wp_parser_delete(m_parser[tid]);
    }
    ::delete[] m_parser;
    ::delete[] m_var;
#endif
}

#endif
