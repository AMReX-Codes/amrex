//BL_COPYRIGHT_NOTICE

//
// $Id: tDir.cpp,v 1.2 1998-10-07 21:18:28 vince Exp $
//

#include <Utility.H>

int
main (int argc, char** argv)
{
    if (argc == 2)
    {
        if (!Utility::UtilCreateDirectory(argv[1], 0755))
        {
            std::cout << "Utility::UtilCreateDirectory() failed!!!\n";
        }
    }
}
