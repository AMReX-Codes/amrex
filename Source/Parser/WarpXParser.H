/* Copyright 2019 Weiqun Zhang
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
#ifndef WARPX_PARSER_H_
#define WARPX_PARSER_H_

#include <array>
#include <vector>
#include <string>
#include <set>

#include <AMReX_REAL.H>

#include "wp_parser_c.h"
#include "wp_parser_y.h"

#ifdef _OPENMP
#include <omp.h>
#endif

template <int N> class GpuParser;

class WarpXParser
{
public:
    WarpXParser (std::string const& func_body);
    WarpXParser () = default;
    ~WarpXParser ();
    void define (std::string const& func_body);

    void setConstant (std::string const& name, amrex::Real c);

    //
    // Option 1: Register every variable to an address provided.
    //           Assign values to external variables.
    //           Call eval().
    void registerVariable (std::string const& name, amrex::Real& var);
    //
    inline amrex::Real eval () const noexcept;

    //
    // Option 2: Register all variables at once. Parser will create
    //               variables internally.
    //           Call eval(...) with variable values.
    void registerVariables (std::vector<std::string> const& names);
    //
    template <typename T, typename... Ts> inline
    amrex::Real eval (T x, Ts... yz) const noexcept;

    void print () const;

    std::string const& expr () const;

    std::set<std::string> symbols () const;

    template <int N> friend class GpuParser;

private:
    void clear ();

    template <typename T> inline
    void unpack (amrex::Real* p, T x) const noexcept;

    template <typename T, typename... Ts> inline
    void unpack (amrex::Real* p, T x, Ts... yz) const noexcept;

    std::string m_expression;
#ifdef _OPENMP
    std::vector<struct wp_parser*> m_parser;
    mutable std::vector<std::array<amrex::Real,16> > m_variables;
    mutable std::vector<std::vector<std::string> > m_varnames;
#else
    struct wp_parser* m_parser = nullptr;
    mutable std::array<amrex::Real,16> m_variables;
    mutable std::vector<std::string> m_varnames;
#endif
};

inline
amrex::Real
WarpXParser::eval () const noexcept
{
#ifdef _OPENMP
    return wp_ast_eval(m_parser[omp_get_thread_num()]->ast,nullptr);
#else
    return wp_ast_eval(m_parser->ast,nullptr);
#endif
}

template <typename T, typename... Ts>
inline
amrex::Real
WarpXParser::eval (T x, Ts... yz) const noexcept
{
#ifdef _OPENMP
    unpack(m_variables[omp_get_thread_num()].data(), x, yz...);
#else
    unpack(m_variables.data(), x, yz...);
#endif
    return eval();
}

template <typename T>
inline
void
WarpXParser::unpack (amrex::Real* p, T x) const noexcept
{
    *p = x;
}

template <typename T, typename... Ts>
inline
void
WarpXParser::unpack (amrex::Real* p, T x, Ts... yz) const noexcept
{
    *p++ = x;
    unpack(p, yz...);
}

#endif
