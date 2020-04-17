#ifndef AMREX_PARTICLEMPIUTIL_H_
#define AMREX_PARTICLEMPIUTIL_H_

#include <map>

#include <AMReX_Vector.H>

namespace amrex {

#ifdef AMREX_USE_MPI    

    Long CountSnds(const std::map<int, Vector<char> >& not_ours, Vector<Long>& Snds);

    Long doHandShake(const std::map<int, Vector<char> >& not_ours,
                     Vector<Long>& Snds, Vector<Long>& Rcvs);

    Long doHandShakeLocal(const std::map<int, Vector<char> >& not_ours,
                          const Vector<int>& neighbor_procs, Vector<Long>& Snds, Vector<Long>& Rcvs);

#endif // AMREX_USE_MPI

}

#endif // include guard
