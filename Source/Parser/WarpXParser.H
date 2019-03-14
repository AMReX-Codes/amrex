#ifndef WARPX_PARSER_H_
#define WARPX_PARSER_H_

#include <vector>
#include <string>
#include <cstdarg>

#include "wp_parser_c.h"
#include "wp_parser_y.h"

#ifdef _OPENMP
#include <omp.h>
#endif


class WarpXParser
{
private:
    std::string m_expression;
    std::vector<struct wp_parser*> m_parser;
    mutable std::vector<std::vector<double> > m_variables;

public:
    WarpXParser (std::string const& func_body);
    WarpXParser () = default;
    ~WarpXParser ();
    void define (std::string const& func_body);

    void setConstant (std::string const& name, double c);

    //
    // Option 1: Register every variable to an address provided.
    //           Assign values to external variables.
    //           Call eval().
    void registerVariable (std::string const& name, double& var);
    //
    inline double eval () const noexcept {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        return wp_ast_eval(m_parser[tid]->ast);
    }

    //
    // Option 2: Register all variables at once. Parser will create
    //               variables internally.
    //           Call eval(...) with variable values.
    void registerVariables (std::vector<std::string> const& names);
    //
    inline double eval (double x, ...) const noexcept {
#ifdef _OPENMP
        const int tid = omp_get_thread_num();
#else
        const int tid = 0;
#endif
        auto& v = m_variables[tid];
        va_list args;
        va_start(args, x);
        v[0] = x;
        for (int i = 1; i < v.size(); ++i) {
            v[i] = va_arg(args, double);
        }
        va_end(args);

        return wp_ast_eval(m_parser[tid]->ast);
    }

    void print () const;

    std::string const& expr () const;

private:
    void clear ();
};

#endif
