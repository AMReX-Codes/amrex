
#include <AMReX_BoxArray.H>
#include <AMReX_Print.H>

using namespace amrex;

extern "C" {

    void amrex_fi_new_boxarray (BoxArray*& ba, int lo[3], int hi[3])
    {
	IntVect small(lo), big(hi);
	ba = new BoxArray(Box(small,big));
    }

    void amrex_fi_new_boxarray_from_bxfarr (BoxArray*& ba, const int* bxs, 
                                            const int nsides, const int ndims, const int nbxs)
    {
        AMREX_ASSERT(nsides == 2);
        AMREX_ASSERT(ndims >= AMREX_SPACEDIM);
        BoxList bl;
        for (int i = 0; i < nbxs; ++i)
        {
            IntVect lo{AMREX_D_DECL( *bxs   , *(bxs+2), *(bxs+4))};
            IntVect hi{AMREX_D_DECL(*(bxs+1), *(bxs+3), *(bxs+5))};
            bl.push_back(Box{lo,hi});
            bxs += 2*ndims;
        }
	ba = new BoxArray(bl);
    }

    void amrex_fi_delete_boxarray (BoxArray* ba)
    {
	delete ba;
    }

    void amrex_fi_clone_boxarray (BoxArray*& bao, const BoxArray* bai)
    {
	delete bao;
	bao = new BoxArray(*bai);
    }

    void amrex_fi_boxarray_maxsize (BoxArray* ba, int sz[])
    {
        IntVect iv(D_DECL(sz[0],sz[1],sz[2]));
        ba->maxSize(iv);
    }

    void amrex_fi_boxarray_get_box (const BoxArray* ba, int i, int lo[3], int hi[3])
    {
        const Box& bx = (*ba)[i];
        const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

    void amrex_fi_print_boxarray (const BoxArray* ba)
    {
	AllPrint() << *ba;
    }
}
