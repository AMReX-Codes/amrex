//BL_COPYRIGHT_NOTICE

//
// $Id: tVisMF.cpp,v 1.2 1997-11-10 21:17:11 lijewski Exp $
//

#include <VisMF.H>

int
main ()
{
    Box bx[] =
    {
        Box(IntVect(D_DECL(0,0,0)), IntVect(D_DECL( 2, 2, 2))),
        Box(IntVect(D_DECL(3,3,3)), IntVect(D_DECL( 8, 8, 8))),
        Box(IntVect(D_DECL(9,9,9)), IntVect(D_DECL(16,16,16)))
    };

    const int NBX = sizeof(bx) / sizeof(bx[0]);

    BoxArray ba(NBX);

    for (int i = 0; i < NBX; i++)
        ba.set(i,bx[i]);

    MultiFab mf(ba, 2, 0);

    for (int i = 0; i < NBX; i++)
        mf[i].setVal(i);

    static const aString mfName = "Spam-n-Eggs";

    VisMF::Write(mf, mfName, VisMF::OneFilePerFab);

    VisMF vmf(mfName);

    for (int i = 0; i < NBX; i++)
    {
        cout << "Fab:\n" << vmf[i] << '\n';
    }
}
