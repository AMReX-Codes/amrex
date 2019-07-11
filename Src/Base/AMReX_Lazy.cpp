#include <AMReX_Lazy.H>

namespace amrex {

namespace Lazy
{
    FuncQue reduction_queue;

    void QueueReduction (Func f)
    {
#ifdef BL_USE_MPI
	reduction_queue.push_back(f);
	const int max_queue_size = 64;
	if (reduction_queue.size() >= max_queue_size)
	    EvalReduction();
#else
	f();
#endif
    }

    void EvalReduction ()
    {
#ifdef BL_USE_MPI
	static int count = 0;
	++count;
	if (count == 1) {
	    for (auto&& f : reduction_queue)
		f();
	    reduction_queue.clear();
            count = 0;
        }
#endif
    }

    void Finalize ()
    {
	EvalReduction();
    }
}

}
