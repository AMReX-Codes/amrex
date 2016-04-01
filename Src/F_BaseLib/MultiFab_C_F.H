#ifndef _MultiFab_C_F_H_
#define _MultiFab_C_F_H_

#include <string>
#include <MultiFab.H>
#include <Geometry.H>

extern "C" {
    void build_layout_from_c (int nb, int dm, const int* lo, const int* hi, 
			      const int* pd_lo, const int* pd_hi,
			      const int* pmask, const int* pmap);
    void destroy_multifab_c ();
    void share_multifab_with_f (const char* fmf_name, int nc, int ng, const int* nodal);
    void share_fab_with_f (int i, const Real* p);
}

class MultiFab_C_to_F
{
public:

    MultiFab_C_to_F (const Geometry& geom,
		     const DistributionMapping& dmap,
		     const BoxArray& ba);

    ~MultiFab_C_to_F ();

    void share (MultiFab& cmf, const std::string& fmf_name);

private:

    static int count;

    struct p_t
    {
	Real* rp;
    };

};


#endif
