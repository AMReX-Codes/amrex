//BL_COPYRIGHT_NOTICE

//
// $Id: Laplacian.cpp,v 1.4 1998-07-29 19:09:54 lijewski Exp $
//

#include <Laplacian.H>

#include <LP_F.H>

void
Laplacian::Fsmooth (MultiFab&       solnL,
                    const MultiFab& rhsL,
                    int             level,
                    int             redBlackFlag)
{
    const BoxArray& bxa = gbox[level];
    OrientationIter oitr;
    const FabSet& f0 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet& f1 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet& f2 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet& f3 = (*undrrelxr[level])[oitr()]; oitr++;
#if (BL_SPACEDIM > 2)
    const FabSet& f4 = (*undrrelxr[level])[oitr()]; oitr++;
    const FabSet& f5 = (*undrrelxr[level])[oitr()]; oitr++;
#endif
    int nc = rhsL.nComp();

    for (MultiFabIterator solnLmfi(solnL); solnLmfi.isValid(); ++solnLmfi)
    {
        DependentMultiFabIterator rhsLmfi(solnLmfi, rhsL);
        DependentFabSetIterator f0fsi(solnLmfi, f0);
        DependentFabSetIterator f1fsi(solnLmfi, f1);
        DependentFabSetIterator f2fsi(solnLmfi, f2);
        DependentFabSetIterator f3fsi(solnLmfi, f3);
#if (BL_SPACEDIM > 2)
        DependentFabSetIterator f4fsi(solnLmfi, f4);
        DependentFabSetIterator f5fsi(solnLmfi, f5);
#endif
        oitr.rewind();
        int gn = solnLmfi.index();
        const Mask& m0 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m1 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m2 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m3 = *maskvals[level][gn][oitr()]; oitr++;
#if (BL_SPACEDIM > 2 )
        const Mask& m4 = *maskvals[level][gn][oitr()]; oitr++;
        const Mask& m5 = *maskvals[level][gn][oitr()]; oitr++;
#endif
        assert(bxa[solnLmfi.index()] == solnLmfi.validbox());

#if (BL_SPACEDIM == 2)
        FORT_GSRB(
            solnLmfi().dataPtr(), 
            ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
            rhsLmfi().dataPtr(), 
            ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
            f0fsi().dataPtr(), 
            ARLIM(f0fsi().loVect()), ARLIM(f0fsi().hiVect()),
            m0.dataPtr(), 
            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
            f1fsi().dataPtr(), 
            ARLIM(f1fsi().loVect()), ARLIM(f1fsi().hiVect()),
            m1.dataPtr(), 
            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
            f2fsi().dataPtr(), 
            ARLIM(f2fsi().loVect()), ARLIM(f2fsi().hiVect()),
            m2.dataPtr(), 
            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
            f3fsi().dataPtr(), 
            ARLIM(f3fsi().loVect()), ARLIM(f3fsi().hiVect()),
            m3.dataPtr(), 
            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
            solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(), &nc,
            h[level], &redBlackFlag);
#endif

#if (BL_SPACEDIM == 3)
        FORT_GSRB(
            solnLmfi().dataPtr(), 
            ARLIM(solnLmfi().loVect()),ARLIM(solnLmfi().hiVect()),
            rhsLmfi().dataPtr(), 
            ARLIM(rhsLmfi().loVect()), ARLIM(rhsLmfi().hiVect()),
            f0fsi().dataPtr(), 
            ARLIM(f0fsi().loVect()), ARLIM(f0fsi().hiVect()),
            m0.dataPtr(), 
            ARLIM(m0.loVect()), ARLIM(m0.hiVect()),
            f1fsi().dataPtr(), 
            ARLIM(f1fsi().loVect()), ARLIM(f1fsi().hiVect()),
            m1.dataPtr(), 
            ARLIM(m1.loVect()), ARLIM(m1.hiVect()),
            f2fsi().dataPtr(), 
            ARLIM(f2fsi().loVect()), ARLIM(f2fsi().hiVect()),
            m2.dataPtr(), 
            ARLIM(m2.loVect()), ARLIM(m2.hiVect()),
            f3fsi().dataPtr(), 
            ARLIM(f3fsi().loVect()), ARLIM(f3fsi().hiVect()),
            m3.dataPtr(), 
            ARLIM(m3.loVect()), ARLIM(m3.hiVect()),
            f4fsi().dataPtr(), 
            ARLIM(f4fsi().loVect()), ARLIM(f4fsi().hiVect()),
            m4.dataPtr(), 
            ARLIM(m4.loVect()), ARLIM(m4.hiVect()),
            f5fsi().dataPtr(), 
            ARLIM(f5fsi().loVect()), ARLIM(f5fsi().hiVect()),
            m5.dataPtr(), 
            ARLIM(m5.loVect()), ARLIM(m5.hiVect()),
            solnLmfi.validbox().loVect(), solnLmfi.validbox().hiVect(), &nc,
            h[level], &redBlackFlag);
#endif
    }
}

void
Laplacian::Fapply (MultiFab&       y,
                   const MultiFab& x,
                   int             level)
{
    int nc = y.nComp();

    for (MultiFabIterator ymfi(y); ymfi.isValid(); ++ymfi)
    {
        DependentMultiFabIterator xmfi(ymfi,  x);

        FORT_ADOTX(ymfi().dataPtr(), 
                   ARLIM(ymfi().loVect()), ARLIM(ymfi().hiVect()),
                   xmfi().dataPtr(), 
                   ARLIM(xmfi().loVect()), ARLIM(xmfi().hiVect()),
                   ymfi.validbox().loVect(), ymfi.validbox().hiVect(), &nc,
                   h[level]);
    }
}

