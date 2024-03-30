#include <AMReX_IParser.H>
#ifdef _WIN32
#define YY_NO_UNISTD_H
#endif
#include <amrex_iparser.lex.h>
#include <amrex_iparser.tab.h>

#include <algorithm>

namespace amrex {

IParser::IParser (std::string const& func_body)
{
    define(func_body);
}

void
IParser::define (std::string const& func_body)
{
    m_data = std::make_shared<Data>();

    if (!func_body.empty()) {
        m_data->m_expression = func_body;
        m_data->m_expression.erase(std::remove(m_data->m_expression.begin(),
                                               m_data->m_expression.end(),'\n'),
                                   m_data->m_expression.end());
        std::string f = m_data->m_expression + "\n";

        YY_BUFFER_STATE buffer = amrex_iparser_scan_string(f.c_str());
        try {
            amrex_iparserparse();
        } catch (const std::runtime_error& e) {
            throw std::runtime_error(std::string(e.what()) + " in IParser expression \""
                                     + m_data->m_expression + "\"");
        }
        m_data->m_iparser = amrex_iparser_new();
        amrex_iparser_delete_buffer(buffer);
    }
}

IParser::Data::~Data ()
{
    m_expression.clear();
    if (m_iparser) { amrex_iparser_delete(m_iparser); }
    if (m_host_executor) { The_Pinned_Arena()->free(m_host_executor); }
#ifdef AMREX_USE_GPU
    if (m_device_executor) { The_Arena()->free(m_device_executor); }
#endif
}

IParser::operator bool () const
{
    return m_data && m_data->m_iparser;
}

void
IParser::setConstant (std::string const& name, int c)
{
    if (m_data && m_data->m_iparser) {
        iparser_setconst(m_data->m_iparser, name.c_str(), c);
    }
}

void
IParser::registerVariables (Vector<std::string> const& vars)
{
    if (m_data && m_data->m_iparser) {
        m_data->m_nvars = static_cast<int>(vars.size());
        for (int i = 0; i < m_data->m_nvars; ++i) {
            iparser_regvar(m_data->m_iparser, vars[i].c_str(), i);
        }
    }
}

void
IParser::print () const
{
    if (m_data && m_data->m_iparser) {
        iparser_print(m_data->m_iparser);
    }
}

int
IParser::depth () const
{
    if (m_data && m_data->m_iparser) {
        return iparser_depth(m_data->m_iparser);
    } else {
        return 0;
    }
}

int
IParser::maxStackSize () const
{
    if (m_data && m_data->m_iparser) {
        return m_data->m_max_stack_size;
    } else {
        return 0;
    }
}

std::string
IParser::expr () const
{
    if (m_data && m_data->m_iparser) {
        return m_data->m_expression;
    } else {
        return std::string{};
    }
}

std::set<std::string>
IParser::symbols () const
{
    if (m_data && m_data->m_iparser) {
        return iparser_get_symbols(m_data->m_iparser);
    } else {
        return std::set<std::string>{};
    }
}

}
