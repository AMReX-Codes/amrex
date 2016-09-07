
#include <MemProfiler.H>
#include <MemProfiler_f.H>

extern "C"
{
    extern long memprof_fab_numdoubles ();
    extern long memprof_fab_numdoubles_hwm ();

    extern long memprof_boxarray_bytes ();
    extern long memprof_boxarray_bytes_hwm ();

    extern long memprof_tilearray_bytes ();
    extern long memprof_tilearray_bytes_hwm ();

    extern long memprof_boxhash_bytes ();
    extern long memprof_boxhash_bytes_hwm ();

    extern long memprof_boxassoc_bytes ();
    extern long memprof_boxassoc_bytes_hwm ();

    extern long memprof_fgassoc_bytes ();
    extern long memprof_fgassoc_bytes_hwm ();

    extern long memprof_syncassoc_bytes ();
    extern long memprof_syncassoc_bytes_hwm ();

    extern long memprof_copyassoc_bytes ();
    extern long memprof_copyassoc_bytes_hwm ();

    extern long memprof_fluxassoc_bytes ();
    extern long memprof_fluxassoc_bytes_hwm ();
}

void 
MemProfiler_f::initialize ()
{
    static bool initialized = false;
    static long sizeof_double = long(sizeof(double));
    if (!initialized) {
	initialized = true;
	MemProfiler::add("F_fab", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_fab_numdoubles()*sizeof_double,
				     memprof_fab_numdoubles_hwm()*sizeof_double};
			 }));
	MemProfiler::add("F_boxarray", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_boxarray_bytes(),
				     memprof_boxarray_bytes_hwm()};
			 }));
	MemProfiler::add("F_tilearray", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_tilearray_bytes(),
				     memprof_tilearray_bytes_hwm()};
			 }));
	MemProfiler::add("F_boxhash", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_boxhash_bytes(),
				     memprof_boxhash_bytes_hwm()};
			 }));
	MemProfiler::add("F_boxassoc", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_boxassoc_bytes(),
				     memprof_boxassoc_bytes_hwm()};
			 }));
	MemProfiler::add("F_fgassoc", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_fgassoc_bytes(),
				     memprof_fgassoc_bytes_hwm()};
			 }));
	MemProfiler::add("F_syncassoc", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_syncassoc_bytes(),
				     memprof_syncassoc_bytes_hwm()};
			 }));
	MemProfiler::add("F_copyassoc", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_copyassoc_bytes(),
				     memprof_copyassoc_bytes_hwm()};
			 }));
	MemProfiler::add("F_fluxassoc", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {memprof_fluxassoc_bytes(),
				     memprof_fluxassoc_bytes_hwm()};
			 }));
    }
}
