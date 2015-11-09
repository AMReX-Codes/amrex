// TODO: cehcking the status of allocate

#include <UArena.H>
#include <ParallelDescriptor.H>

void*
UArena::alloc (std::size_t _sz)
{
    if (ParallelDescriptor::TeamSize() > 1) {
	return upcxx::allocate(_sz);
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
