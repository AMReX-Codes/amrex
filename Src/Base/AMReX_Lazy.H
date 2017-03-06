#ifndef BL_LAZY_H
#define BL_LAZY_H

#include <vector>
#include <functional>
#include <algorithm>

namespace amrex {
namespace Lazy
{
    typedef typename std::function<void()> Func;
    typedef typename std::vector<Func> FuncQue;

    extern FuncQue reduction_queue;

    void QueueReduction (Func);
    void EvalReduction ();
    
    void Finalize ();
}
}

#endif
