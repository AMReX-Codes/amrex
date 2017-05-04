#include <cstring>
#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_REAL.H>
#include <AMReX_Print.H>

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

    int amrex_parmparse_get_counts (ParmParse* pp, const char* name)
    {
	return pp->countval(name);
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

    void amrex_parmparse_get_string (ParmParse* pp, const char* name, char* v, int len)
    {
        std::string b;
        pp->get(name, b);
        std::strncpy(v, b.c_str(), len);
    }

    void amrex_parmparse_get_intarr (ParmParse* pp, const char* name, int v[], int len)
    {
	Array<int> r;
	pp->getarr(name, r);
	for (int i = 0; i < len; ++i) {
	    v[i] = r[i];
	}
    }

    void amrex_parmparse_get_realarr (ParmParse* pp, const char* name, Real v[], int len)
    {
	Array<Real> r;
	pp->getarr(name, r);
	for (int i = 0; i < len; ++i) {
	    v[i] = r[i];
	}
    }

    void amrex_parmparse_get_stringarr (ParmParse* pp, const char* name, char* v, int len, int n)
    {
        std::vector<std::string> b;
        pp->getarr(name, b);
        BL_ASSERT(n == b.size());
        for (int i = 0; i < n; ++i) {
            std::strncpy(v, b[i].c_str(), len);
            v += len;
        }
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

    void amrex_parmparse_query_string (ParmParse* pp, const char* name, char* v, int len)
    {
      std::string b;
      pp->query(name, b);
      std::strncpy(v, b.c_str(), len);
    }
}
