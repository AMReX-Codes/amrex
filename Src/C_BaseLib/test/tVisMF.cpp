//BL_COPYRIGHT_NOTICE

//
// $Id: tVisMF.cpp,v 1.4 1997-11-11 01:18:16 lijewski Exp $
//

#include <stdlib.h>

#include <VisMF.H>

static
void
Write_N_Read (const MultiFab& mf,
              const aString&  mf_name,
              VisMF::How      how)
{
    switch (how)
    {
    case VisMF::OneFilePerCPU:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerCPU); break;
    case VisMF::OneFilePerFab:
        VisMF::Write(mf, mf_name, VisMF::OneFilePerFab); break;
    default:
        BoxLib::Error("Bad case in switch");
    }

    VisMF vmf(mf_name);

    for (int i = 0, N = vmf.length(); i < N; i++)
    {
        cout << "Fab:\n" << vmf[i] << '\n';
    }

    aString cmd("/bin/rm -f ");

    cmd += mf_name;
    cmd += "*";

    system(cmd.c_str());
}

int
main ()
{
    StartParallel(2);

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

    MultiFab mf(ba, 2, 1);

    for (int i = 0; i < NBX; i++)
    {
        mf[i].setVal(i+1);
    }

    static const aString mf_name = "Spam-n-Eggs";

    Write_N_Read (mf, mf_name, VisMF::OneFilePerCPU);

    Write_N_Read (mf, mf_name, VisMF::OneFilePerFab);

    EndParallel();
}
