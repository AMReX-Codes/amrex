#ifndef BL_CCSE_MPP
#define BL_CCSE_MPP

//
// $Id: ccse-mpi.H,v 1.11 2001-09-25 22:45:41 lijewski Exp $
//

//
// This header encapsulates the MPI header.
//
// By placing <mpi.h> here we make it easier to recompile with or w/o MPI.
//

#ifdef BL_USE_MPI
#include <winstd.H>

#include <mpi.h>

//
// Template function to return MPI_Datatype depending on argument type.
//
// We only implement this function for standard MPI data types.
//
template<class T> MPI_Datatype mpi_data_type(T*);
//
// And the appropriate specializations for the basic MPI types.
//
template<> inline MPI_Datatype mpi_data_type(char*)           { return MPI_CHAR;   }
template<> inline MPI_Datatype mpi_data_type(short*)          { return MPI_SHORT;  }
template<> inline MPI_Datatype mpi_data_type(int*)            { return MPI_INT;    }
template<> inline MPI_Datatype mpi_data_type(long*)           { return MPI_LONG;   }
template<> inline MPI_Datatype mpi_data_type(float*)          { return MPI_FLOAT;  }
template<> inline MPI_Datatype mpi_data_type(double*)         { return MPI_DOUBLE; }
template<> inline MPI_Datatype mpi_data_type(unsigned char*)  { return MPI_UNSIGNED_CHAR;   }
template<> inline MPI_Datatype mpi_data_type(unsigned short*) { return MPI_UNSIGNED_SHORT;  }
template<> inline MPI_Datatype mpi_data_type(unsigned int*)   { return MPI_UNSIGNED;    }
template<> inline MPI_Datatype mpi_data_type(unsigned long*)  { return MPI_UNSIGNED_LONG;   }

namespace ParallelDescriptor
{
    template <class T>
    struct Mpi_typemap
    {
	static MPI_Datatype type();
    };
    template <> MPI_Datatype Mpi_typemap<char>::type();
    template <> MPI_Datatype Mpi_typemap<short>::type();
    template <> MPI_Datatype Mpi_typemap<int>::type();
    template <> MPI_Datatype Mpi_typemap<long>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned char>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned short>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned int>::type();
    template <> MPI_Datatype Mpi_typemap<unsigned long>::type();
    template <> MPI_Datatype Mpi_typemap<float>::type();
    template <> MPI_Datatype Mpi_typemap<double>::type();

    struct Max         { static MPI_Op op(); };
    struct Min         { static MPI_Op op(); };
    struct Sum         { static MPI_Op op(); };
    struct Prod        { static MPI_Op op(); };
    struct Logical_And { static MPI_Op op(); };
    struct Boolean_And { static MPI_Op op(); };
    struct Logical_Or  { static MPI_Op op(); };
    struct Boolean_Or  { static MPI_Op op(); };
    struct Logical_XOr { static MPI_Op op(); };
    struct Boolean_XOr { static MPI_Op op(); };
}


#else

typedef int MPI_Op;
typedef int MPI_Comm;
typedef int MPI_Datatype;
struct MPI_Status {};
typedef int MPI_Request;
const int MPI_REQUEST_NULL = 0;

#endif /*BL_USE_MPI*/

#endif /*BL_CCSE_MPP*/
