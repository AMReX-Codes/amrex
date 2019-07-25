#include <GpuParser.H>

GpuParser::GpuParser (WarpXParser const& wp)
{
#ifdef AMREX_USE_GPU

    struct wp_parser* a_wp = wp.m_parser;
    m_gpu_parser.sz_mempool = wp_ast_size((struct wp_node*)a_wp);
    m_gpu_parser.p_root = (struct wp_node*)
        amrex::The_Managed_Arena()->alloc(m_gpu_parser.sz_mempool);
    m_gpu_parser.p_free = m_gpu_parser.p_root;
    // 0: don't free the source
    m_gpu_parser.ast = wp_parser_ast_dup(&m_gpu_parser, a_wp->ast, 0);
    wp_parser_regvar_gpu(&m_gpu_parser, "x", 0);
    wp_parser_regvar_gpu(&m_gpu_parser, "y", 1);
    wp_parser_regvar_gpu(&m_gpu_parser, "z", 2);

    m_cpu_parser.sz_mempool = wp_ast_size((struct wp_node*)a_wp);
    m_cpu_parser.p_root = (struct wp_node*)
        amrex::The_Managed_Arena()->alloc(m_cpu_parser.sz_mempool);
    m_cpu_parser.p_free = m_cpu_parser.p_root;
    // 0: don't free the source
    m_cpu_parser.ast = wp_parser_ast_dup(&m_cpu_parser, a_wp->ast, 0);
    wp_parser_regvar(&m_cpu_parser, "x", &(m_var.x));
    wp_parser_regvar(&m_cpu_parser, "y", &(m_var.y));
    wp_parser_regvar(&m_cpu_parser, "z", &(m_var.z));
    
#else

#ifdef _OPENMP
    nthreads = omp_get_max_threads();
#else
    nthreads = 1;
#endif

    m_parser = ::new struct wp_parser*[nthreads];
    m_var = ::new amrex::XDim3[nthreads];

    for (int tid = 0; tid < nthreads; ++tid)
    {
#ifdef _OPENMP
        m_parser[tid] = wp_parser_dup(wp.m_parser[tid]);
#else
        m_parser[tid] = wp_parser_dup(wp.m_parser);
#endif
        wp_parser_regvar(m_parser[tid], "x", &(m_var[tid].x));
        wp_parser_regvar(m_parser[tid], "y", &(m_var[tid].y));
        wp_parser_regvar(m_parser[tid], "z", &(m_var[tid].z));
    }

#endif
}

void
GpuParser::clear ()
{
#ifdef AMREX_USE_GPU
    amrex::The_Managed_Arena()->free(m_gpu_parser.ast);
    amrex::The_Managed_Arena()->free(m_cpu_parser.ast);
#else
    for (int tid = 0; tid < nthreads; ++tid)
    {
        wp_parser_delete(m_parser[tid]);
    }
    ::delete[] m_parser;
    ::delete[] m_var;
#endif
}

