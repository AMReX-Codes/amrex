#ifndef _CACHE_H_
#define _CACHE_H_

#  include "interface.H"
#  include "boundary.H"

#ifdef BL_FORT_USE_UNDERSCORE
#  define FORT_FFCCPY1     fccpy1_
#  define FORT_FFCCPY2     fccpy2_
#else
#  define FORT_FFCCPY1     FCCPY1
#  define FORT_FFCCPY2     FCCPY2
#endif

extern "C"
{
#if (BL_SPACEDIM == 2)
    void FORT_FFCCPY1(Real*, const Real*, const int*, int*, int*, int*, int*, int*);
#else
    void FORT_FFCCPY2(Real*, const Real*, const int*, int*, int*, int*, int*, int*, int*, int*, int*);
#endif
}

class copy_cache 
{
public:
    copy_cache(int Nsets, Real *Dptr, Real *Sptr);
    
    // sync cache
    copy_cache(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy);
    
    // border cache (only w == 1 case currently supported)
    copy_cache(MultiFab& r, const level_interface& lev_interface, const amr_boundary_class* bdy, int w, bool hg_terrain);
    
    ~copy_cache(); 
#if (BL_SPACEDIM == 2)
    void set(int i, int Dstart, int Sstart, int Dstrid, int Sstrid, int Nvals) 
    {
	if (i < 0 || i >= nsets)
	    BoxLib::Error("copy_cache::set---out of range");
	dstart[i] = Dstart;
	sstart[i] = Sstart;
	dstrid[i] = Dstrid;
	sstrid[i] = Sstrid;
	nvals[i]  = Nvals;
    }
#else
    void set(int i, int Dstart, int Sstart, int Dstrid1, int Dstrid2, int Sstrid1, int Sstrid2, int Nvals1, int Nvals2) 
    {
	if (i < 0 || i >= nsets)
	    BoxLib::Error("copy_cache::set---out of range");
	dstart[i] = Dstart;
	sstart[i] = Sstart;
	if (Nvals1 >= Nvals2) 
	{
	    dstrid1[i] = Dstrid1;
	    dstrid2[i] = Dstrid2;
	    sstrid1[i] = Sstrid1;
	    sstrid2[i] = Sstrid2;
	    nvals1[i]  = Nvals1;
	    nvals2[i]  = Nvals2;
	}
	else 
	{
	    dstrid1[i] = Dstrid2;
	    dstrid2[i] = Dstrid1;
	    sstrid1[i] = Sstrid2;
	    sstrid2[i] = Sstrid1;
	    nvals1[i]  = Nvals2;
	    nvals2[i]  = Nvals1;
	}
    }
#endif
    
    void run() const ;

    const Real* dataPtr() const;
    Real* dataPtr();

//private:
    // data are public so that users may set up their own calls
    int nsets;
    Real *dptr, *sptr;
    int *dstart, *sstart;
#if (BL_SPACEDIM == 2)
    int *dstrid, *sstrid;
    int *nvals;
#else
    int *dstrid1, *dstrid2, *sstrid1, *sstrid2;
    int *nvals1, *nvals2;
#endif
};

inline const Real* copy_cache::dataPtr() const
{
    return dptr;
}


inline Real* copy_cache::dataPtr()
{
    return dptr;
}

inline copy_cache::~copy_cache() 
{
    delete [] dstart;
}


inline void copy_cache::run() const 
{
#if (BL_SPACEDIM == 2)
    FORT_FFCCPY1(dptr, sptr, &nsets, dstart, sstart, dstrid, sstrid, nvals);
#else
    FORT_FFCCPY2(dptr, sptr, &nsets, dstart, sstart, dstrid1, dstrid2, sstrid1, sstrid2, nvals1, nvals2);
#endif
}

class unroll_cache 
{
public:
    unroll_cache(int Nsets, Real *Ptr)
	: nsets(Nsets), ptr(Ptr) 
    {
#if (BL_SPACEDIM == 2)
	start = new int[3 * nsets];
	strid = start + nsets;
	nvals = start + 2 * nsets;
#else
	start  = new int[4 * nsets];
	strid1 = start + nsets;
	strid2 = start + 2 * nsets;
	nvals  = start + 3 * nsets;
#endif
    }
    unroll_cache(MultiFab&);
    ~unroll_cache() 
    {
	delete [] start;
    }
    
#if (BL_SPACEDIM == 2)
    void set(int i, int Start, int Strid, int Nvals) 
    {
	if (i < 0 || i >= nsets)
	    BoxLib::Error("unroll_cache::set---out of range");
	start[i] = Start;
	strid[i] = Strid;
	nvals[i] = Nvals;
    }
#else
    void set(int i, int Start, int Strid1, int Strid2, int Nvals) 
    {
	if (i < 0 || i >= nsets)
	    BoxLib::Error("unroll_cache::set---out of range");
	start[i]  = Start;
	strid1[i] = Strid1;
	strid2[i] = Strid2;
	nvals[i]  = Nvals;
    }
#endif
//private:
    // data are public so that users may set up their own calls
    const int nsets;
    Real* ptr;
    int *start;
#if (BL_SPACEDIM == 2)
    int *strid;
#else
    int *strid1, *strid2;
#endif
    int *nvals;
    
};

#endif

