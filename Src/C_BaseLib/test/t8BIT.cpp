
//
// $Id: t8BIT.cpp,v 1.2 2000-10-02 20:52:40 lijewski Exp $
//
// A simple program to read in a MultiFab and write out in 8BIT format.
//

#include <VisMF.H>

int
main (int argc, char** argv)
{
    argc--; argv++;

    FArrayBox::setFormat(FABio::FAB_8BIT);

    for (int i = 0; i < argc; i++)
    {
        cout << "Transforming " << argv[i] << " ... " << flush;

        aString name = argv[i];

        MultiFab mf;

        VisMF::Read(mf, name);

        VisMF::Write(mf, name, VisMF::OneFilePerCPU, true);

        cout << "done" << endl;
    }
}
