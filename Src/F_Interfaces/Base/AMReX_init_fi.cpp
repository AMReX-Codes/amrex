#include <AMReX.H>
#include <AMReX_Array.H>
#include <AMReX_Vector.H>

#include <string>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <cstring>
#include <iostream>

extern "C"
{
    using amrex_void_cfun = void (*)();

    void amrex_fi_init (char const* cmd, int fcomm, int arg_parmparse, amrex_void_cfun proc_parmparse)
    {
        std::istringstream is(cmd);
        amrex::Vector<std::string> argv_string(std::istream_iterator<std::string>{is},
                                              std::istream_iterator<std::string>{  });

        auto argc = static_cast<int>(argv_string.size());
        auto** argv = (char**)malloc(argc*sizeof(char*));
        for (int i = 0; i < argc; ++i)
        {
            argv[i] = (char*)malloc(argv_string[i].size()+1);
            std::strncpy(argv[i], argv_string[i].c_str(), argv_string[i].size()+1);
        }

#ifdef BL_USE_MPI
        MPI_Comm ccomm = MPI_Comm_f2c(fcomm);
#else
        amrex::ignore_unused(fcomm);
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
