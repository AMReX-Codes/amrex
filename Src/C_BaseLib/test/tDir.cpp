
//
// $Id: tDir.cpp,v 1.4 2001-07-22 18:25:48 car Exp $
//

#include <Utility.H>

int
main (int argc, char** argv)
{
    if (argc == 2)
    {
        if (!BoxLib::UtilCreateDirectory(argv[1], 0755))
        {
            std::cout << "Utility::UtilCreateDirectory() failed!!!\n";
        }
    }
}
