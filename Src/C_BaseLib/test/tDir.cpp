
//
// $Id: tDir.cpp,v 1.3 2000-10-02 20:52:40 lijewski Exp $
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
