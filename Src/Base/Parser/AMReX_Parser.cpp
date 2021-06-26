#include <AMReX_Parser.H>
#ifdef _WIN32
#define YY_NO_UNISTD_H
#endif
#include <amrex_parser.lex.h>
#include <amrex_parser.tab.h>

#include <algorithm>

namespace amrex {

Parser::Parser (std::string const& func_body)
{
    define(func_body);
}

void
Parser::define (std::string const& func_body)
{
    clear();

    m_expression = func_body;
    m_expression.erase(std::remove(m_expression.begin(),m_expression.end(),'\n'),
                       m_expression.end());
    std::string f = m_expression + "\n";

    YY_BUFFER_STATE buffer = amrex_parser_scan_string(f.c_str());
    amrex_parserparse();
    m_parser = amrex_parser_new();
    amrex_parser_delete_buffer(buffer);
}

Parser::~Parser ()
{
    clear();
}

void
Parser::clear ()
{
    m_expression.clear();
    if (m_parser) { amrex_parser_delete(m_parser); }
    m_parser = nullptr;
    if (m_host_executor) { The_Pinned_Arena()->free(m_host_executor); }
    m_host_executor = nullptr;
#ifdef AMREX_USE_GPU
    if (m_device_executor) { The_Arena()->free(m_device_executor); }
    m_device_executor = nullptr;
#endif
}

void
Parser::setConstant (std::string const& name, amrex::Real c)
{
    parser_setconst(m_parser, name.c_str(), c);
}

void
Parser::registerVariables (Vector<std::string> const& vars)
{
    m_nvars = vars.size();
    for (int i = 0; i < m_nvars; ++i) {
        parser_regvar(m_parser, vars[i].c_str(), i);
    }
}

void
Parser::print () const
{
    parser_print(m_parser);
}

int
Parser::depth () const
{
    return parser_depth(m_parser);
}

int
Parser::maxStackSize () const
{
    return m_max_stack_size;
}

std::string const&
Parser::expr () const
{
    return m_expression;
}

std::set<std::string>
Parser::symbols () const
{
    return parser_get_symbols(m_parser);
}

}
