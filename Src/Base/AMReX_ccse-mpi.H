
#ifndef BL_CCSE_MPI_H_
#define BL_CCSE_MPI_H_
#include <AMReX_Config.H>

//
// This header encapsulates the MPI header.
//
// By placing <mpi.h> here we make it easier to recompile with or w/o MPI.
//

#ifdef BL_USE_MPI

#include <mpi.h>

namespace amrex::ParallelDescriptor
{
    using lull_t = unsigned long long[8];

    template <class T>
    //! Communication datatype (note: this structure also works without MPI)
    struct Mpi_typemap
    {
        static MPI_Datatype type();
    };
    //
    // Template functions to return MPI_Datatype depending on argument type.
    //
    template <> MPI_Datatype Mpi_typemap<char>::type();
    template <> MPI_Datatype Mpi_typemap<short>::type();
    template <> MPI_Datatype Mpi_typemap<int>::type();
    template <> MPI_Datatype Mpi_typemap<long>::type();
    template <> MPI_Datatype Mpi_typemap<long long>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned char>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned short>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned int>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned long>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned long long>::type();
    template <> MPI_Datatype Mpi_typemap<float>::type();
    template <> MPI_Datatype Mpi_typemap<double>::type();
    template <> MPI_Datatype Mpi_typemap<lull_t>::type();
}

#else

using MPI_Op       = int;
using MPI_Comm     = int;
using MPI_Group    = int;
using MPI_Datatype = int;
using MPI_Request  = int;
struct MPI_Status {};
static constexpr int MPI_DATATYPE_NULL = 0;
static constexpr int MPI_REQUEST_NULL  = 0;
static constexpr int MPI_COMM_WORLD    = 0;
static constexpr int MPI_COMM_NULL     = 0;
static constexpr int MPI_GROUP_NULL    = 0;
static constexpr int MPI_PROC_NULL     = 0;
static constexpr int MPI_MAX_PROCESSOR_NAME = 64;

static constexpr int MPI_MINLOC = 0;
static constexpr int MPI_MAXLOC = 0;

namespace amrex::ParallelDescriptor
{
    template <class T>
    //! Communication datatype (note: this structure also works without MPI)
    struct Mpi_typemap
    {
        static MPI_Datatype type();
    };
}

#endif /*BL_USE_MPI*/

#endif
