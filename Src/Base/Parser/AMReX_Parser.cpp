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
    m_data = std::make_shared<Data>();

    if (!func_body.empty()) {
        m_data->m_expression = func_body;
        m_data->m_expression.erase(std::remove(m_data->m_expression.begin(),
                                               m_data->m_expression.end(),'\n'),
                                   m_data->m_expression.end());
        std::string f = m_data->m_expression + "\n";

        YY_BUFFER_STATE buffer = amrex_parser_scan_string(f.c_str());
        try {
            amrex_parserparse();
        } catch (const std::runtime_error& e) {
            amrex_parser_delete_buffer(buffer); // delete buffer allocated by bison
            amrex_parser_delete_ptrs();         // delete ptrs allocated by amrex
            throw std::runtime_error(std::string(e.what()) + " in Parser expression \""
                                     + m_data->m_expression + "\"");
        }
        m_data->m_parser = amrex_parser_new();
        amrex_parser_delete_buffer(buffer);
        m_ufs = parser_get_user_functions(m_data->m_parser);
    }
}

Parser::Data::~Data ()
{
    m_expression.clear();
    if (m_parser) { amrex_parser_delete(m_parser); }
    if (m_host_executor) {
        if (m_use_arena) {
            The_Pinned_Arena()->free(m_host_executor);
        } else {
            std::free(m_host_executor);
        }
    }
#ifdef AMREX_USE_GPU
    if (m_device_executor) { The_Arena()->free(m_device_executor); }
#endif
}

Parser::operator bool () const
{
    return m_data && m_data->m_parser;
}

void
Parser::setConstant (std::string const& name, double c)
{
    if (m_data && m_data->m_parser) {
        parser_setconst(m_data->m_parser, name.c_str(), c);
    }
}

void
Parser::registerVariables (Vector<std::string> const& vars)
{
    m_vars = vars;
    if (m_data && m_data->m_parser) {
        m_data->m_nvars = static_cast<int>(vars.size());
        for (int i = 0; i < m_data->m_nvars; ++i) {
            parser_regvar(m_data->m_parser, vars[i].c_str(), i);
        }
    }
}

void
Parser::registerUserFn1 (std::string const& name, ParserUserFn1 fh, ParserUserFn1 fd)
{
    register_user_fn<1>(name,fh,fd);
}

void
Parser::registerUserFn2 (std::string const& name, ParserUserFn2 fh, ParserUserFn2 fd)
{
    register_user_fn<2>(name,fh,fd);
}

void
Parser::registerUserFn3 (std::string const& name, ParserUserFn3 fh, ParserUserFn3 fd)
{
    register_user_fn<3>(name,fh,fd);
}

void
Parser::registerUserFn4 (std::string const& name, ParserUserFn4 fh, ParserUserFn4 fd)
{
    register_user_fn<4>(name,fh,fd);
}

void
Parser::print () const
{
    if (m_data && m_data->m_parser) {
        parser_print(m_data->m_parser);
    }
}

int
Parser::depth () const
{
    if (m_data && m_data->m_parser) {
        return parser_depth(m_data->m_parser);
    } else {
        return 0;
    }
}

int
Parser::maxStackSize () const
{
    if (m_data && m_data->m_parser) {
        return m_data->m_max_stack_size;
    } else {
        return 0;
    }
}

std::string
Parser::expr () const
{
    if (m_data && m_data->m_parser) {
        return m_data->m_expression;
    } else {
        return std::string{};
    }
}

std::set<std::string>
Parser::symbols () const
{
    if (m_data && m_data->m_parser) {
        return parser_get_symbols(m_data->m_parser);
    } else {
        return std::set<std::string>{};
    }
}

std::map<std::string,int> const&
Parser::userFunctions () const
{
    return m_ufs;
}

void
Parser::printExe () const
{
    if (m_data->m_host_executor) {
        parser_exe_print(m_data->m_host_executor, m_vars, m_data->m_locals);
    }
}

}
