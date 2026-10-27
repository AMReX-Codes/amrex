#include <AMReX_BArena.H>
#include <AMReX.H>

void*
amrex::BArena::alloc (std::size_t sz_)
{
    sz_ = Arena::align(sz_);
    auto sys_alloc = [] (std::size_t sz) {
#ifndef _WIN32
        return std::aligned_alloc(align_size, sz);
#else
        return _aligned_malloc(sz, align_size);
#endif
    };
    void* pt = sys_alloc(sz_);
    if (pt == nullptr && sz_ > 0) {
        // the pooled host arenas may hold enough unused memory
        auto release = [] (Arena* a) -> std::size_t {
            return (a->isHostAccessible() && !a->isManaged()) ? a->freeUnused() : 0;
        };
        std::size_t nfreed = release(The_Comms_Arena());
        if (The_Pinned_Arena() != The_Comms_Arena()) {
            nfreed += release(The_Pinned_Arena());
        }
        if (nfreed > 0) { pt = sys_alloc(sz_); }
    }
    if (pt == nullptr && sz_ > 0) {
        amrex::ErrorStream() << "Out of CPU memory: got nullptr from std::aligned_alloc\n";
        throw_out_of_memory("CPU memory", sz_, "std::aligned_alloc returned nullptr");
    }
    m_profiler.profile_alloc(pt, sz_);
    return pt;
}

void
amrex::BArena::free (void* pt)
{
    m_profiler.profile_free(pt);
#ifndef _WIN32
    std::free(pt);
#else
    _aligned_free(pt);
#endif
}

bool
amrex::BArena::isDeviceAccessible () const
{
    return false;
}

bool
amrex::BArena::isHostAccessible () const
{
    return true;
}

bool
amrex::BArena::isManaged () const
{
    return false;
}

bool
amrex::BArena::isDevice () const
{
    return false;
}

bool
amrex::BArena::isPinned () const
{
    return false;
}
