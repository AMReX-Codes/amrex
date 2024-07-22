#ifndef BL_LAZY_H
#define BL_LAZY_H
#include <AMReX_Config.H>

#include <vector>
#include <functional>
#include <algorithm>

namespace amrex {
namespace Lazy
{
    using Func = std::function<void()>;
    using FuncQue = std::vector<Func>;

    extern FuncQue reduction_queue;

    void QueueReduction (Func);
    void EvalReduction ();

    void Finalize ();
}
}

#endif
