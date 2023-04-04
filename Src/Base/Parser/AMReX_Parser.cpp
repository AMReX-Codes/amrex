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
            throw std::runtime_error(std::string(e.what()) + " in Parser expression \""
                                     + m_data->m_expression + "\"");
        }
        m_data->m_parser = amrex_parser_new();
        amrex_parser_delete_buffer(buffer);
    }
}

Parser::Data::~Data ()
{
    m_expression.clear();
    if (m_parser) { amrex_parser_delete(m_parser); }
    if (m_host_executor) { The_Pinned_Arena()->free(m_host_executor); }
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
    if (m_data && m_data->m_parser) {
        m_data->m_nvars = static_cast<int>(vars.size());
        for (int i = 0; i < m_data->m_nvars; ++i) {
            parser_regvar(m_data->m_parser, vars[i].c_str(), i);
        }
    }
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

}
