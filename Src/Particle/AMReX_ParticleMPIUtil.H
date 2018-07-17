#ifndef AMREX_PARTICLEMPIUTIL_H_
#define AMREX_PARTICLEMPIUTIL_H_

#include <map>

#include <AMReX_Vector.H>

namespace amrex {

#ifdef BL_USE_MPI    

    long CountSnds(const std::map<int, Vector<char> >& not_ours, Vector<long>& Snds);

    long doHandShake(const std::map<int, Vector<char> >& not_ours,
                     Vector<long>& Snds, Vector<long>& Rcvs);

    long doHandShakeLocal(const std::map<int, Vector<char> >& not_ours,
                          const Vector<int>& neighbor_procs, Vector<long>& Snds, Vector<long>& Rcvs);

#endif // BL_USE_MPI

}

#endif // include guard
