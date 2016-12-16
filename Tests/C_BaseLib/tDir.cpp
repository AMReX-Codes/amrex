
#include <AMReX_Utility.H>

int
main (int argc, char** argv)
{
    if (argc == 2)
    {
        if (!amrex::UtilCreateDirectory(argv[1], 0755))
        {
            std::cout << "Utility::UtilCreateDirectory() failed!!!\n";
        }
    }
}
