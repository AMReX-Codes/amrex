#ifndef BL_CCSE_MPP
#define BL_CCSE_MPP

//
// $Id: ccse-mpi.H,v 1.7 2001-07-17 23:02:30 lijewski Exp $
//

//
// This header encapsulates the MPI header.
//
// By placing <mpi.h> here we make it easier to recompile with or w/o MPI.
//

#ifdef BL_USE_MPI
#include <mpi.h>
#ifdef WIN32
#undef small
#undef min
#undef max
#endif

#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif
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
#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif
#else
#ifdef BL_NAMESPACE
namespace BL_NAMESPACE
{
#endif
typedef int MPI_Op;
typedef int MPI_Comm;
#ifdef BL_NAMESPACE
}
#ifndef BL_NO_USING_DIRECTIVE
using namespace BL_NAMESPACE;
#endif
#endif
#endif /*BL_USE_MPI*/

#endif /*BL_CCSE_MPP*/
