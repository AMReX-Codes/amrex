
#include "WarpXParser.H"

WarpXParser::WarpXParser (std::string const& func_body)
{
    define(func_body);
}

void
WarpXParser::define (std::string const& func_body)
{
    clear();

    m_expression = func_body;

#ifdef _OPENMP
    int nthreads = omp_get_max_threads();
#else
    int nthreads = 1;
#endif
    m_parser.resize(nthreads);

    std::string f = m_expression + "\n";
    m_parser[0] = wp_c_parser_new(f.c_str());

#ifdef _OPENMP
#pragma omp parallel
    {
        int tid = omp_get_thread_num();
        if (tid > 0) {
            m_parser[tid] = wp_parser_dup(m_parser[0]);
        }
    }
#endif
}

WarpXParser::~WarpXParser ()
{
    clear();
}

void
WarpXParser::clear ()
{
    if (!m_parser.empty())
    {
#ifdef _OPENMP
#pragma omp parallel
        {
            int tid = omp_get_thread_num();
            wp_parser_delete(m_parser[tid]);
        }
#else
        wp_parser_delete(m_parser[0]);
#endif
    }
    m_expression.clear();
    m_parser.clear();
    m_variables.clear();
}

void
WarpXParser::registerVariable (std::string const& name, double& var)
{
    // We assume this is called inside OMP parallel region
#ifdef _OPENMP
    const int tid = omp_get_thread_num();
#else
    const int tid = 0;
#endif
    wp_parser_regvar(m_parser[tid], name.c_str(), &var);
}

// This must be called outside OpenMP parallel region.
void
WarpXParser::registerVariables (std::vector<std::string> const& names)
{
#ifdef _OPENMP
    const int nthreads = omp_get_max_threads();
#else
    const int nthreads = 1;
#endif
    m_variables.resize(nthreads);
    const int nnames = names.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < nthreads; ++i) {
        m_variables[i].resize(nnames);
        for (int j = 0; j < nnames; ++i) {
            wp_parser_regvar(m_parser[i], names[j].c_str(), &(m_variables[i][j]));
        }
    }
}

void
WarpXParser::setConstant (std::string const& name, double c)
{
    // We don't know if this is inside OMP parallel region or not
#ifdef _OPENMP
    bool in_parallel = omp_in_parallel();
#pragma omp parallel if (!in_parallel)
    {
        wp_parser_setconst(m_parser[omp_get_thread_num()], name.c_str(), c);
    }
#else
    wp_parser_setconst(m_parser[0], name.c_str(), c);
#endif
}

void
WarpXParser::print () const
{
#ifdef _OPENMP
#pragma omp critical(warpx_parser_pint)
    wp_ast_print(m_parser[omp_get_thread_num()]->ast);
#else
    wp_ast_print(m_parser[0]->ast);
#endif
}

std::string const&
WarpXParser::expr () const
{
    return m_expression;
}
