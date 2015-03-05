#ifndef BL_LAZY_H
#define BL_LAZY_H

#include <deque>
#include <functional>
#include <algorithm>

namespace Lazy
{
    typedef std::function<void()> Func;
    typedef std::vector<Func> FuncQue;

    extern FuncQue reduction_queue;

    void QueueReduction (Func);
    void EvalReduction ();
    
    void Finalize ();
}

#endif
