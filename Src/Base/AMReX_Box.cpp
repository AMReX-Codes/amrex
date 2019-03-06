
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
AllGatherBoxes (Vector<Box>& bxs)
{
#ifdef BL_USE_MPI
    // cell centered boxes only!
    const auto szof_bx = Box::linearSize();

    const long count = bxs.size() * static_cast<long>(szof_bx);
    const auto& countvec = ParallelDescriptor::Gather(count, ParallelDescriptor::IOProcessorNumber());
    
    long count_tot = 0L;
    Vector<long> offset(countvec.size(),0L);
    if (ParallelDescriptor::IOProcessor())
    {
        count_tot = countvec[0];
        for (int i = 1, N = offset.size(); i < N; ++i) {
            offset[i] = offset[i-1] + countvec[i-1];
            count_tot += countvec[i];
        }
    }

    ParallelDescriptor::Bcast(&count_tot, 1, ParallelDescriptor::IOProcessorNumber());

    if (count_tot == 0) return;

    Vector<char> send_buffer(count);
    char* psend = (count > 0) ? send_buffer.data() : nullptr;
    char* p = psend;
    for (const auto& b : bxs) {
        b.linearOut(p);
        p += szof_bx;
    }

    Vector<char> recv_buffer(count_tot);
    ParallelDescriptor::Gatherv(psend, count, recv_buffer.data(), countvec, offset, ParallelDescriptor::IOProcessorNumber());

    ParallelDescriptor::Bcast(recv_buffer.data(), count_tot, ParallelDescriptor::IOProcessorNumber());

    const long nboxes_tot = count_tot/szof_bx;
    bxs.resize(nboxes_tot);

    p = recv_buffer.data();
    for (auto& b : bxs) {
        b.linearIn(p);
        p += szof_bx;
    }
#endif
}

}
