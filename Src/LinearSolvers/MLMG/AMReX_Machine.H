#ifndef AMREX_MACHINE_H
#define AMREX_MACHINE_H

namespace amrex {
namespace machine {

void Initialize (); // called in amrex::Initialize()

// find the best topologically close neighborhood of ranks
// returns a vector of global or local rank IDs based on flag_local_ranks
Vector<int> find_best_nbh (int rank_n, bool flag_local_ranks = false);

}}

#endif
