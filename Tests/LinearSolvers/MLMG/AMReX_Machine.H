#ifndef AMREX_MACHINE_H
#define AMREX_MACHINE_H

namespace amrex {
namespace machine {

void init (); // must be called by ALL ranks in the job
Vector<int> find_best_nbh (int rank_n);

}}

#endif
