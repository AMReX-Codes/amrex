
#include <winstd.H>

#include <cstring>
#include <BaseFab.H>
#include <BArena.H>
#include <CArena.H>
#include <Thread.H>
#include <SPECIALIZE_F.H>

long BoxLib::total_bytes_allocated_in_fabs = 0;

long BoxLib::total_bytes_allocated_in_fabs_hwm = 0;

int BoxLib::BF_init::m_cnt = 0;

static ThreadSpecificData<Arena>* arena = 0;

BoxLib::BF_init::BF_init ()
{
    if (m_cnt++ == 0)
        arena = new ThreadSpecificData<Arena>;
}

BoxLib::BF_init::~BF_init ()
{
    if (--m_cnt == 0)
        delete arena;
}

Arena*
BoxLib::ResetArena (Arena* newarena)
{
    BL_ASSERT(newarena != 0);

    Arena* oldarena = arena->get();

    BL_ASSERT(oldarena != 0);

    arena->set(newarena);

    return oldarena;
}

Arena*
BoxLib::The_Arena ()
{
    BL_ASSERT(arena != 0);

    Arena* a = arena->get();

    if (a == 0)
    {
#if defined(BL_THREADS) || defined(BL_COALESCE_FABS)
        arena->set(a = new CArena);
#else
        arena->set(a = new BArena);
#endif
    }

    return a;
}

#if !( defined(WIN32))
template<>
void
BaseFab<Real>::performCopy (const BaseFab<Real>& src,
                            const Box&           srcbox,
                            int                  srccomp,
                            const Box&           destbox,
                            int                  destcomp,
                            int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    if (destbox == domain && srcbox == src.box())
    {
        Real*       data_dst = dataPtr(destcomp);
        const Real* data_src = src.dataPtr(srccomp);

        for (long i = 0, N = numcomp*numpts; i < N; i++)
        {
            *data_dst++ = *data_src++;
        }
    }
    else
    {
        const int* destboxlo  = destbox.loVect();
        const int* destboxhi  = destbox.hiVect();
        const int* _th_plo    = loVect();
        const int* _th_phi    = hiVect();
        const int* _x_lo      = srcbox.loVect();
        const int* _x_hi      = srcbox.hiVect(); 
        const int* _x_plo     = src.loVect();
        const int* _x_phi     = src.hiVect();
        Real*       _th_p     = dataPtr(destcomp);
        const Real* _x_p      = src.dataPtr(srccomp);

        FORT_FASTCOPY(_th_p,
                      ARLIM(_th_plo),
                      ARLIM(_th_phi),
                      D_DECL(destboxlo[0],destboxlo[1],destboxlo[2]),
                      D_DECL(destboxhi[0],destboxhi[1],destboxhi[2]),
                      _x_p,
                      ARLIM(_x_plo),
                      ARLIM(_x_phi),
                      D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                      &numcomp);
    }
}

template<>
void
BaseFab<Real>::performSetVal (Real       val,
                              const Box& bx,
                              int        comp,
                              int        ncomp)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

    Real* data = dataPtr(comp);

    if (bx == domain)
    {
        for (long i = 0, N = ncomp*numpts; i < N; i++)
        {
            *data++ = val;
        }
    }
    else
    {
        const int* _box_lo = bx.loVect();            
        const int* _box_hi = bx.hiVect();            
        const int* _th_plo = loVect();                           
        const int* _th_phi = hiVect();

        FORT_FASTSETVAL(&val,
                        _box_lo,
                        _box_hi,
                        data,
                        ARLIM(_th_plo),
                        ARLIM(_th_phi),
                        &ncomp);
    }
}

template<>
Real
BaseFab<Real>::norm (const Box& bx,
                     int        p,
                     int        comp,
                     int        ncomp) const
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(comp >= 0 && comp + ncomp <= nvar);

    const int* _box_lo = bx.loVect();            
    const int* _box_hi = bx.hiVect();            
    const int* _datalo = loVect();                           
    const int* _datahi = hiVect();

    const Real* _data = dataPtr(comp);

    Real nrm = 0;

    if (p == 0)
    {
        FORT_FASTZERONORM(_data,
                          ARLIM(_datalo),
                          ARLIM(_datahi),
                          _box_lo,
                          _box_hi,
                          &ncomp,
                          &nrm);
    }
    else if (p == 1)
    {
        FORT_FASTONENORM(_data,
                         ARLIM(_datalo),
                         ARLIM(_datahi),
                         _box_lo,
                         _box_hi,
                         &ncomp,
                         &nrm);
    }
    else
    {
        BoxLib::Error("BaseFab<Real>::norm(): only p == 0 or p == 1 are supported");
    }

    return nrm;
}

template<>
BaseFab<Real>&
BaseFab<Real>::plus (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const int* destboxlo  = destbox.loVect();
    const int* destboxhi  = destbox.hiVect();
    const int* _th_plo    = loVect();
    const int* _th_phi    = hiVect();
    const int* _x_lo      = srcbox.loVect();
    const int* _x_hi      = srcbox.hiVect(); 
    const int* _x_plo     = src.loVect();
    const int* _x_phi     = src.hiVect();
    Real*       _th_p     = dataPtr(destcomp);
    const Real* _x_p      = src.dataPtr(srccomp);

    FORT_FASTPLUS(_th_p,
                  ARLIM(_th_plo),
                  ARLIM(_th_phi),
                  D_DECL(destboxlo[0],destboxlo[1],destboxlo[2]),
                  D_DECL(destboxhi[0],destboxhi[1],destboxhi[2]),
                  _x_p,
                  ARLIM(_x_plo),
                  ARLIM(_x_phi),
                  D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                  &numcomp);

    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::mult (const BaseFab<Real>& src,
                     const Box&           srcbox,
                     const Box&           destbox,
                     int                  srccomp,
                     int                  destcomp,
                     int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const int* destboxlo  = destbox.loVect();
    const int* destboxhi  = destbox.hiVect();
    const int* _th_plo    = loVect();
    const int* _th_phi    = hiVect();
    const int* _x_lo      = srcbox.loVect();
    const int* _x_hi      = srcbox.hiVect(); 
    const int* _x_plo     = src.loVect();
    const int* _x_phi     = src.hiVect();
    Real*       _th_p     = dataPtr(destcomp);
    const Real* _x_p      = src.dataPtr(srccomp);

    FORT_FASTMULT(_th_p,
                  ARLIM(_th_plo),
                  ARLIM(_th_phi),
                  D_DECL(destboxlo[0],destboxlo[1],destboxlo[2]),
                  D_DECL(destboxhi[0],destboxhi[1],destboxhi[2]),
                  _x_p,
                  ARLIM(_x_plo),
                  ARLIM(_x_phi),
                  D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                  &numcomp);

    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::minus (const BaseFab<Real>& src,
                      const Box&           srcbox,
                      const Box&           destbox,
                      int                  srccomp,
                      int                  destcomp,
                      int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const int* destboxlo  = destbox.loVect();
    const int* destboxhi  = destbox.hiVect();
    const int* _th_plo    = loVect();
    const int* _th_phi    = hiVect();
    const int* _x_lo      = srcbox.loVect();
    const int* _x_hi      = srcbox.hiVect(); 
    const int* _x_plo     = src.loVect();
    const int* _x_phi     = src.hiVect();
    Real*       _th_p     = dataPtr(destcomp);
    const Real* _x_p      = src.dataPtr(srccomp);

    FORT_FASTMINUS(_th_p,
                   ARLIM(_th_plo),
                   ARLIM(_th_phi),
                   D_DECL(destboxlo[0],destboxlo[1],destboxlo[2]),
                   D_DECL(destboxhi[0],destboxhi[1],destboxhi[2]),
                   _x_p,
                   ARLIM(_x_plo),
                   ARLIM(_x_phi),
                   D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                   &numcomp);

    return *this;
}

template<>
BaseFab<Real>&
BaseFab<Real>::divide (const BaseFab<Real>& src,
                       const Box&           srcbox,
                       const Box&           destbox,
                       int                  srccomp,
                       int                  destcomp,
                       int                  numcomp)
{
    BL_ASSERT(destbox.ok());
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());

    const int* destboxlo  = destbox.loVect();
    const int* destboxhi  = destbox.hiVect();
    const int* _th_plo    = loVect();
    const int* _th_phi    = hiVect();
    const int* _x_lo      = srcbox.loVect();
    const int* _x_hi      = srcbox.hiVect(); 
    const int* _x_plo     = src.loVect();
    const int* _x_phi     = src.hiVect();
    Real*       _th_p     = dataPtr(destcomp);
    const Real* _x_p      = src.dataPtr(srccomp);

    FORT_FASTDIVIDE(_th_p,
                    ARLIM(_th_plo),
                    ARLIM(_th_phi),
                    D_DECL(destboxlo[0],destboxlo[1],destboxlo[2]),
                    D_DECL(destboxhi[0],destboxhi[1],destboxhi[2]),
                    _x_p,
                    ARLIM(_x_plo),
                    ARLIM(_x_phi),
                    D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                    &numcomp);

    return *this;
}
#endif
