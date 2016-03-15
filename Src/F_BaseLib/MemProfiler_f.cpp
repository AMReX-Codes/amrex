
#include <MemProfiler.H>
#include <MemProfiler_f.H>

extern "C"
{
    extern long memprof_fab_numdoubles ();
    extern long memprof_fab_numdoubles_hwm ();
}

void 
MemProfiler_f::initialize ()
{
    static bool initialized = false;
    if (!initialized) {
	initialized = true;
	MemProfiler::add("Fab_f", [] () -> MemProfiler::MemInfo {
		return {memprof_fab_numdoubles()*sizeof(double),
			memprof_fab_numdoubles_hwm()*sizeof(double)};
	    });
    }
}
