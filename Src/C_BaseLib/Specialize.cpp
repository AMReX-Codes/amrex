//
// $Id: Specialize.cpp,v 1.2 1999-05-10 17:18:47 car Exp $
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
    BLassert(src.box().contains(srcbox));
    BLassert(box().contains(destbox));
    BLassert(destbox.sameSize(srcbox));
    BLassert(srccomp >= 0 && srccomp+numcomp <= src.nComp());
    BLassert(destcomp >= 0 && destcomp+numcomp <= nComp());
 
    Box _subbox_(box()); 
    _subbox_ &= destbox; 
    BLassert((srcbox).sameSize(_subbox_)); 
 
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
    BLassert(domain.contains(bx));
    BLassert(ns >= 0 && ns + num <= nvar);

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
