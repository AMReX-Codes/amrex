
#ifndef BL_CCSE_MPI_H_
#define BL_CCSE_MPI_H_
#include <AMReX_Config.H>

/**
 * \file AMReX_ccse-mpi.H
 * \brief MPI datatype helpers and fallbacks for builds without MPI.
 *
 * Including this header centralizes all direct references to <mpi.h>.
 */

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

    /**
     * \brief Helper that maps a C++ type to its corresponding MPI_Datatype.
     *
     * Specializations are provided for all fundamental types AMReX sends.
     * This wrapper is available even in non-MPI builds (where MPI_Datatype
     * is an int placeholder), so templated code can compile unconditionally.
     */
    template <class T>
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

namespace amrex::mpidatatypes {}

#else

namespace amrex::mpidatatypes {

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

}

namespace amrex {
    using namespace mpidatatypes;
}

namespace amrex::ParallelDescriptor
{
    /**
     * \brief Helper that maps a C++ type to its corresponding MPI_Datatype.
     *
     * In non-MPI builds the returned value is a dummy placeholder so code
     * that relies on the interface can still compile.
     */
    template <class T>
    struct Mpi_typemap
    {
        static MPI_Datatype type();
    };
}

#endif /*BL_USE_MPI*/

#endif
