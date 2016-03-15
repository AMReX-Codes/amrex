
#include <MemProfiler.H>
#include <MemProfiler_f.H>

extern "C"
{
    extern long memprof_fab_numdoubles ();
    extern long memprof_fab_numdoubles_hwm ();

    extern long memprof_boxarray_bytes ();
    extern long memprof_boxarray_bytes_hwm ();

    extern long memprof_boxassoc_bytes ();
    extern long memprof_boxassoc_bytes_hwm ();

    extern long memprof_copyassoc_bytes ();
    extern long memprof_copyassoc_bytes_hwm ();
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
	MemProfiler::add("BoxArray_f", [] () -> MemProfiler::MemInfo {
		return {memprof_boxarray_bytes(),
			memprof_boxarray_bytes_hwm()};
	    });
	MemProfiler::add("Boxassoc_f", [] () -> MemProfiler::MemInfo {
		return {memprof_boxassoc_bytes(),
			memprof_boxassoc_bytes_hwm()};
	    });
	MemProfiler::add("Copyassoc_f", [] () -> MemProfiler::MemInfo {
		return {memprof_copyassoc_bytes(),
			memprof_copyassoc_bytes_hwm()};
	    });
    }
}
