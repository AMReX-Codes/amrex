#ifndef AMREX_REDUCE_OP_TYPE_H_
#define AMREX_REDUCE_OP_TYPE_H_

namespace amrex {

/// \ingroup amrex_collectives
enum class ReduceOpType : int {
    max = 0,
    min,
    sum,
    lor,
    land
};

}

#endif
