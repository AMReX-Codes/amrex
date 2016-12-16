#include <cstring>
#include <AMReX_ParmParse.H>

extern "C"
{
    void fi_new_parmparse (ParmParse*& pp, const char* name)
    {
	pp = new ParmParse(std::string(name));
    }

    void fi_delete_parmparse (ParmParse* pp)
    {
	delete pp;
    }

    void fi_parmparse_get_int (ParmParse* pp, const char* name, int* v)
    {
	pp->get(name, *v);
    }

    void fi_parmparse_get_double (ParmParse* pp, const char* name, double* v)
    {
	pp->get(name, *v);
    }

    void fi_parmparse_get_bool (ParmParse* pp, const char* name, int* v)
    {
	bool b;
	pp->get(name, b);
	*v = b;
    }

    void fi_parmparse_get_string (ParmParse* pp, const char* name, char* v, int* len)
    {
      std::string b;
      pp->get(name, b);
      std::strncpy(v, b.c_str(), *len);
    }

    void fi_parmparse_query_int (ParmParse* pp, const char* name, int* v)
    {
	pp->query(name, *v);
    }

    void fi_parmparse_query_double (ParmParse* pp, const char* name, double* v)
    {
	pp->query(name, *v);
    }

    void fi_parmparse_query_bool (ParmParse* pp, const char* name, int* v)
    {
	bool b = *v;
	pp->query(name, b);
	*v = b;
    }

    void fi_parmparse_query_string (ParmParse* pp, const char* name, char* v, int* len)
    {
      std::string b;
      pp->query(name, b);
      std::strncpy(v, b.c_str(), *len);
    }
}
