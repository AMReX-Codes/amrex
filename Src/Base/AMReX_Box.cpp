
#include <iostream>
#include <limits>

#include <AMReX_BLassert.H>
#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Print.H>
#include <AMReX_ParallelDescriptor.H>

namespace amrex {

//
// I/O functions.
//

std::ostream&
operator<< (std::ostream& os,
            const Box&    b)
{
    os << '('
       << b.smallEnd() << ' '
       << b.bigEnd()   << ' '
       << b.type()
       << ')';

    if (os.fail())
        amrex::Error("operator<<(ostream&,Box&) failed");

    return os;
}

//
// Moved out of Utility.H
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            Box&          b)
{
    IntVect lo, hi, typ;

    is >> std::ws;
    char c;
    is >> c;

    if (c == '(')
    {
        is >> lo >> hi;
	is >> c;
	// Read an optional IndexType
	is.putback(c);
	if ( c == '(' )
	{
	    is >> typ;
	}
        is.ignore(BL_IGNORE_MAX,')');
    }
    else if (c == '<')
    {
	is.putback(c);
        is >> lo >> hi;
	is >> c;
	// Read an optional IndexType
	is.putback(c);
	if ( c == '<' )
	{
	    is >> typ;
	}
        //is.ignore(BL_IGNORE_MAX,'>');
    }
    else
    {
        amrex::Error("operator>>(istream&,Box&): expected \'(\'");
    }

    b = Box(lo,hi,typ);

    if (is.fail())
        amrex::Error("operator>>(istream&,Box&) failed");

    return is;
}

BoxCommHelper::BoxCommHelper (const Box& bx, int* p_)
    : p(p_)
{
    if (p == 0) {
	v.resize(3*AMREX_SPACEDIM);
	p = &v[0];
    }

    AMREX_D_EXPR(p[0]               = bx.smallend[0],
	   p[1]               = bx.smallend[1],
	   p[2]               = bx.smallend[2]);
    AMREX_D_EXPR(p[0+AMREX_SPACEDIM]   = bx.bigend[0],
	   p[1+AMREX_SPACEDIM]   = bx.bigend[1],
	   p[2+AMREX_SPACEDIM]   = bx.bigend[2]);
    const IntVect& typ = bx.btype.ixType();
    AMREX_D_EXPR(p[0+AMREX_SPACEDIM*2] = typ[0],
	   p[1+AMREX_SPACEDIM*2] = typ[1],
	   p[2+AMREX_SPACEDIM*2] = typ[2]);
}

void
AllGatherBoxes (Vector<Box>& bxs, int n_extra_reserve)
{
#ifdef BL_USE_MPI

#if 0
    // In principle, MPI_Allgather/MPI_Allgatherv should not be slower than
    // MPI_Gather/MPI_Gatherv followed by MPI_Bcast.  But that's not true on Summit.
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    const int count = bxs.size();
    Vector<int> countvec(ParallelContext::NProcsSub());
    MPI_Allgather(&count, 1, MPI_INT, countvec.data(), 1, MPI_INT, comm);

    Vector<int> offset(countvec.size(),0);
    Long count_tot = countvec[0];
    for (int i = 1, N = offset.size(); i < N; ++i) {
        offset[i] = offset[i-1] + countvec[i-1];
        count_tot += countvec[i];
    }

    if (count_tot == 0) return;

    if (count_tot > static_cast<Long>(std::numeric_limits<int>::max())) {
        amrex::Abort("AllGatherBoxes: not many boxes");
    }

    Vector<Box> recv_buffer;
    recv_buffer.reserve(count_tot+n_extra_reserve);
    recv_buffer.resize(count_tot);
    MPI_Allgatherv(bxs.data(), count, ParallelDescriptor::Mpi_typemap<Box>::type(),
                   recv_buffer.data(), countvec.data(), offset.data(),
                   ParallelDescriptor::Mpi_typemap<Box>::type(), comm);

    std::swap(bxs,recv_buffer);
#else
    MPI_Comm comm = ParallelContext::CommunicatorSub();
    const int root = ParallelContext::IOProcessorNumberSub();
    const int myproc = ParallelContext::MyProcSub();
    const int nprocs = ParallelContext::NProcsSub();
    const int count = bxs.size();
    Vector<int> countvec(nprocs);
    MPI_Gather(&count, 1, MPI_INT, countvec.data(), 1, MPI_INT, root, comm);

    Long count_tot = 0L;
    Vector<int> offset(countvec.size(),0);
    if (myproc == root) {
        count_tot = countvec[0];
        for (int i = 1, N = offset.size(); i < N; ++i) {
            offset[i] = offset[i-1] + countvec[i-1];
            count_tot += countvec[i];
        }
    }

    MPI_Bcast(&count_tot, 1, MPI_INT, root, comm);

    if (count_tot == 0) return;

    if (count_tot > static_cast<Long>(std::numeric_limits<int>::max())) {
        amrex::Abort("AllGatherBoxes: not many boxes");
    }

    Vector<Box> recv_buffer;
    recv_buffer.reserve(count_tot+n_extra_reserve);
    recv_buffer.resize(count_tot);
    MPI_Gatherv(bxs.data(), count, ParallelDescriptor::Mpi_typemap<Box>::type(),
                recv_buffer.data(), countvec.data(), offset.data(),
                ParallelDescriptor::Mpi_typemap<Box>::type(), root, comm);
    MPI_Bcast(recv_buffer.data(), count_tot, ParallelDescriptor::Mpi_typemap<Box>::type(),
              root, comm);

    std::swap(bxs,recv_buffer);
#endif

#else
    amrex::ignore_unused(bxs,n_extra_reserve);
#endif
}

}
