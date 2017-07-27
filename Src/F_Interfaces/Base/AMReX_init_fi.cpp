#include <AMReX.H>
#include <AMReX_Array.H>

#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cstring>

extern "C"
{
    void amrex_fi_init (char* cmd, int fcomm, int arg_parmparse,
                        amrex::PTR_TO_VOID_FUNC proc_parmparse)
    {
        std::istringstream is(cmd);
        amrex::Array<std::string> argv_string(std::istream_iterator<std::string>{is},
                                              std::istream_iterator<std::string>{  });

        int argc = argv_string.size();
        char** argv = (char**)malloc(argc*sizeof(char*));
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = (char*)malloc(argv_string[i].size()+1);
            strcpy(argv[i], argv_string[i].c_str());
        }

#ifdef BL_USE_MPI
        MPI_Comm ccomm = MPI_Comm_f2c(fcomm);
#else
        int ccomm = 0;
#endif
        amrex::Initialize(argc, argv, arg_parmparse, ccomm, proc_parmparse);

        for (int i = 0; i < argc; ++i)
        {
            free(argv[i]);
        }
        free(argv);
    }

    void amrex_fi_finalize ()
    {
        amrex::Finalize();
    }
}
