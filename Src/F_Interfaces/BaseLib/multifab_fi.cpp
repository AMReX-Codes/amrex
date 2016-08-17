
#include <MultiFab.H>
#include <Geometry.H>

extern "C" {

    void fi_new_multifab (MultiFab*& mf, BoxArray*& bao, const BoxArray* bai, 
			  int nc, int ng, const int* nodal)
    {
	mf = new MultiFab(*bai, nc, ng, Fab_allocate, IntVect(nodal));
	bao = (BoxArray*)&(mf->boxArray());
    }

    void fi_delete_multifab (MultiFab* mf)
    {
	delete mf;
    }

    void fi_multifab_dataptr (MultiFab* mf, MFIter* mfi, double*& dp, int lo[3], int hi[3])
    {
	FArrayBox& fab = (*mf)[*mfi];
	dp = fab.dataPtr();
	const Box& bx = fab.box();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

    double fi_multifab_min(const MultiFab* mf, int comp, int nghost)
    {
	return mf->min(comp,nghost);
    }

    double fi_multifab_max(const MultiFab* mf, int comp, int nghost)
    {
	return mf->max(comp,nghost);
    }

    double fi_multifab_norm0(const MultiFab* mf, int comp)
    {
	return mf->norm0(comp);
    }

    double fi_multifab_norm1(const MultiFab* mf, int comp)
    {
	return mf->norm1(comp);
    }

    double fi_multifab_norm2(const MultiFab* mf, int comp)
    {
	return mf->norm2(comp);
    }

    void fi_multifab_fill_boundary (MultiFab* mf, const Geometry* geom, 
				    int c, int nc, int cross)
    {
	mf->FillBoundary(c, nc, cross, geom.periodicity());
    }

    // MFIter routines

    void fi_new_mfiter (MFIter*& mfi, MultiFab* mf, int tiling)
    {
	mfi = new MFIter(*mf, (bool)tiling);
    }

    void fi_delete_mfiter (MFIter* mfi)
    {
	delete mfi;
    }

    void fi_increment_mfiter (MFIter* mfi, int* isvalid)
    {
	++(*mfi);
	*isvalid = mfi->isValid();
    }

    void fi_mfiter_is_valid (MFIter* mfi, int* isvalid)
    {
	*isvalid = mfi->isValid();
    }

    void fi_mfiter_tilebox (MFIter* mfi, int lo[3], int hi[3])
    {
	const Box& bx = mfi->tilebox();
	const int* lov = bx.loVect();
	const int* hiv = bx.hiVect();
	for (int i = 0; i < BL_SPACEDIM; ++i) {
	    lo[i] = lov[i];
	    hi[i] = hiv[i];
	}
    }

}

