//
// $Id: Specialize.cpp,v 1.3 1999-05-10 18:54:23 car Exp $
//

#ifdef BL_USE_SPECIALIZE

#include <BaseFab.H>
//
// Explicit <Real> specializations for performCopy() and performSetVal().
//
#include <SPECIALIZE_F.H>

#ifdef BL_SPECIALIZE_SYNTAX
template <>
#endif
void
BaseFab<Real>::performCopy (const BaseFab<Real>& src,
                            const Box&           srcbox,
                            int                  srccomp,
                            const Box&           destbox,
                            int                  destcomp,
                            int                  numcomp)
{
    BL_ASSERT(src.box().contains(srcbox));
    BL_ASSERT(box().contains(destbox));
    BL_ASSERT(destbox.sameSize(srcbox));
    BL_ASSERT(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BL_ASSERT(destcomp >= 0 && destcomp+numcomp <= nComp());
 
    Box _subbox_(box()); 
    _subbox_ &= destbox; 
    BL_ASSERT((srcbox).sameSize(_subbox_)); 
 
    if (_subbox_.ok())
    { 
        const int* _subbox_lo = _subbox_.loVect(); // dest subbox to fill
        const int* _subbox_hi = _subbox_.hiVect(); 

        const int* _th_plo    = loVect();          // bounds of dest array
        const int* _th_phi    = hiVect();
 
        const int* _x_lo      = (srcbox).loVect(); // src subbox to fill from
        const int* _x_hi      = (srcbox).hiVect(); 
        const int* _x_plo     = (src).loVect();    // bounds of src array
        const int* _x_phi     = (src).hiVect();

        Real* _th_p      = dataPtr(destcomp);      // dest array
        const Real* _x_p = (src).dataPtr(srccomp); // src array

        FORT_FASTCOPY(_th_p,
                      ARLIM(_th_plo),
                      ARLIM(_th_phi),
                      D_DECL(_subbox_lo[0],_subbox_lo[1],_subbox_lo[2]),
                      D_DECL(_subbox_hi[0],_subbox_hi[1],_subbox_hi[2]),
                      _x_p,
                      ARLIM(_x_plo),
                      ARLIM(_x_phi),
                      D_DECL(_x_lo[0],_x_lo[1],_x_lo[2]),
                      D_DECL(_x_hi[0],_x_hi[1],_x_hi[2]),
                      numcomp);
    }
}

#ifdef BL_SPECIALIZE_SYNTAX
template <>
#endif
void
BaseFab<Real>::performSetVal (Real       val,
                              const Box& bx,
                              int        ns,
                              int        num)
{
    BL_ASSERT(domain.contains(bx));
    BL_ASSERT(ns >= 0 && ns + num <= nvar);

    const int* _box_lo = bx.loVect();            
    const int* _box_hi = bx.hiVect();            

    const int* _th_plo = loVect();                           
    const int* _th_phi = hiVect();

    Real* _th_p = dataPtr(ns);

    FORT_FASTSETVAL(&val,
                    _box_lo,
                    _box_hi,
                    _th_p,
                    ARLIM(_th_plo),
                    ARLIM(_th_phi),
                    num);
}
#endif /*BL_USE_SPECIALIZE*/
