#include <cstring>
#include <AMReX_ParmParse.H>
#include <AMReX_REAL.H>

using namespace amrex;

extern "C"
{
    void amrex_new_parmparse (ParmParse*& pp, const char* name)
    {
	pp = new ParmParse(std::string(name));
    }

    void amrex_delete_parmparse (ParmParse* pp)
    {
	delete pp;
    }

    void amrex_parmparse_get_int (ParmParse* pp, const char* name, int* v)
    {
	pp->get(name, *v);
    }

    void amrex_parmparse_get_real (ParmParse* pp, const char* name, Real* v)
    {
	pp->get(name, *v);
    }

    void amrex_parmparse_get_bool (ParmParse* pp, const char* name, int* v)
    {
	bool b;
	pp->get(name, b);
	*v = b;
    }

    void amrex_parmparse_get_string (ParmParse* pp, const char* name, char* v, int* len)
    {
      std::string b;
      pp->get(name, b);
      std::strncpy(v, b.c_str(), *len);
    }

    void amrex_parmparse_query_int (ParmParse* pp, const char* name, int* v)
    {
	pp->query(name, *v);
    }

    void amrex_parmparse_query_real (ParmParse* pp, const char* name, Real* v)
    {
	pp->query(name, *v);
    }

    void amrex_parmparse_query_bool (ParmParse* pp, const char* name, int* v)
    {
	bool b = *v;
	pp->query(name, b);
	*v = b;
    }

    void amrex_parmparse_query_string (ParmParse* pp, const char* name, char* v, int* len)
    {
      std::string b;
      pp->query(name, b);
      std::strncpy(v, b.c_str(), *len);
    }
}
