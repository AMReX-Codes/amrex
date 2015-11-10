#include <UArena.H>
#include <ParallelDescriptor.H>
#include <BaseFab.H>

void*
UArena::alloc (std::size_t _sz)
{
    if (ParallelDescriptor::TeamSize() > 1) {
	auto p = upcxx::allocate(_sz);
	if (p == nullptr) {
	    std::cout << "===== Proc. " << ParallelDescriptor::MyProc() << " =====\n";
	    std::cout << "   Failed to allocate " << _sz << " bytes global address space memory!\n";
	    std::cout << "   Please try to increase the GASNET_MAX_SEGSIZE environment variable.\n";
	    std::cout << "   For example, export GASNET_MAX_SEGSIZE=512MB\n";
	    std::cout << "   Total Bytes Allocated in Fabs: " << BoxLib::TotalBytesAllocatedInFabs();
	    std::cout << "   Highest Watermark in Fabs: " << BoxLib::TotalBytesAllocatedInFabsHWM();
	    std::cout << std::endl;
	    BoxLib::Abort("UArena: upcxx::allocate failed");
	}
	return p;
    } else {
	return ::operator new(_sz);
    }
}

void
UArena::free (void* pt)
{
    if (ParallelDescriptor::TeamSize() > 1) {
	upcxx::deallocate(pt);
    } else {
	::operator delete(pt);
    }
}
