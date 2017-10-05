#include <cstring>
#include <AMReX_ParmParse.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>
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

    void amrex_parmparse_get_string (ParmParse* pp, const char* name, char*& v, int* len)
    {
        std::string b;
        pp->get(name, b);
        *len = b.size() + 1;
        v = new char[*len];
        std::strncpy(v, b.c_str(), *len);
    }

    void amrex_parmparse_delete_cp_char (char**v, int len)
    {
        for (int i = 0; i < len; ++i) {
            delete[] v[i];
        }
    }

    void amrex_parmparse_get_intarr (ParmParse* pp, const char* name, int v[], int len)
    {
	Vector<int> r;
	pp->getarr(name, r);
	for (int i = 0; i < len; ++i) {
	    v[i] = r[i];
	}
    }

    void amrex_parmparse_get_realarr (ParmParse* pp, const char* name, Real v[], int len)
    {
	Vector<Real> r;
	pp->getarr(name, r);
	for (int i = 0; i < len; ++i) {
	    v[i] = r[i];
	}
    }

    void amrex_parmparse_get_stringarr (ParmParse* pp, const char* name, char** v, int* sv, int n)
    {
        std::vector<std::string> b;
        pp->getarr(name, b);
        BL_ASSERT(n == static_cast<int>(b.size()));
        for (int i = 0; i < n; ++i) {
            sv[i] = b[i].size() + 1;
            v[i] = new char[sv[i]];
            std::strncpy(v[i], b[i].c_str(), sv[i]);
        }
    }

    int amrex_parmparse_query_int (ParmParse* pp, const char* name, int* v)
    {
	return pp->query(name, *v);
    }

    int amrex_parmparse_query_real (ParmParse* pp, const char* name, Real* v)
    {
	return pp->query(name, *v);
    }

    int amrex_parmparse_query_bool (ParmParse* pp, const char* name, int* v)
    {
	bool b;
	if (pp->query(name, b)) {
            *v = b;
            return 1;
        } else {
            return 0;
        }
    }

    int amrex_parmparse_query_string (ParmParse* pp, const char* name, char*& v, int* len)
    {
      std::string b;
      int r = pp->query(name, b);
      *len = b.size() + 1;
      v = new char[*len];
      std::strncpy(v, b.c_str(), *len);
      return r;
    }

    void amrex_parmparse_add_int (ParmParse* pp, const char* name, int v)
    {
        pp->add(name,v);
    }

    void amrex_parmparse_add_real (ParmParse* pp, const char* name, Real v)
    {
        pp->add(name,v);
    }
    
    void amrex_parmparse_add_bool (ParmParse* pp, const char* name, int v)
    {
        pp->add(name,static_cast<bool>(v));
    }

    void amrex_parmparse_add_string (ParmParse* pp, const char* name, const char* v)
    {
        pp->add(name, std::string{v});
    }

    void amrex_parmparse_add_intarr (ParmParse* pp, const char* name, const int v[], int len)
    {
        pp->addarr(name, std::vector<int>(v, v+len));
    }

    void amrex_parmparse_add_realarr (ParmParse* pp, const char* name, const Real v[], int len)
    {
        pp->addarr(name, std::vector<Real>(v, v+len));
    }

    void amrex_parmparse_add_stringarr (ParmParse* pp, const char* name, const char* v, int len)
    {
        std::vector<std::string> vs;
        vs.reserve(len);
        const char* p = v;
        for (int i = 0; i < len; ++i) {
            vs.push_back(p);
            p += vs[i].size()+1;
        }
        pp->addarr(name, vs);
    }
}
