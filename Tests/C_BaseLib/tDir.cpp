//BL_COPYRIGHT_NOTICE

//
// $Id: tDir.cpp,v 1.1 1997-11-12 21:57:28 lijewski Exp $
//

#include <Utility.H>

int
main (int argc, char** argv)
{
    if (argc == 2)
    {
        if (!Utility::CreateDirectory(argv[1], 0755))
        {
            std::cout << "Utility::CreateDirectory() failed!!!\n";
        }
    }
}
