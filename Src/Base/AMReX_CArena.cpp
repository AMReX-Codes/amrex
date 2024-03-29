
#include <AMReX_CArena.H>
#include <AMReX_BLassert.H>
#include <AMReX_Gpu.H>
#include <AMReX_ParallelReduce.H>

#ifdef AMREX_TINY_PROFILING
#include <AMReX_TinyProfiler.H>
#else
namespace amrex {
    struct MemStat {};
}
#endif

#include <utility>
#include <cstring>
#include <iostream>

namespace amrex {

CArena::CArena (std::size_t hunk_size, ArenaInfo info)
    : m_hunk(align(hunk_size == 0 ? DefaultHunkSize : hunk_size))
{
    arena_info = info;
    BL_ASSERT(m_hunk >= hunk_size);
    BL_ASSERT(m_hunk%Arena::align_size == 0);
}

CArena::~CArena ()
{
    for (auto const& a : m_alloc) {
        deallocate_system(a.first, a.second);
    }

#ifdef AMREX_TINY_PROFILING
    if (m_do_profiling) {
        TinyProfiler::DeregisterArena(m_profiling_stats);
    }
#endif
}

void*
CArena::alloc (std::size_t nbytes)
{
    std::lock_guard<std::mutex> lock(carena_mutex);
    nbytes = Arena::align(nbytes == 0 ? 1 : nbytes);
    return alloc_protected(nbytes);
}

void*
CArena::alloc_protected (std::size_t nbytes)
{
    MemStat* stat = nullptr;
#ifdef AMREX_TINY_PROFILING
    if (m_do_profiling) {
        stat = TinyProfiler::memory_alloc(nbytes, m_profiling_stats);
    }
#endif

    if (static_cast<Long>(m_used+nbytes) >= arena_info.release_threshold) {
        freeUnused_protected();
    }

    //
    // Find node in freelist at lowest memory address that'll satisfy request.
    //
    auto free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it) {
        if ((*free_it).size() >= nbytes) {
            break;
        }
    }

    void* vp = nullptr;

    if (free_it == m_freelist.end())
    {
        const std::size_t N = nbytes < m_hunk ? m_hunk : nbytes;

        vp = allocate_system(N);

        m_used += N;

        m_alloc.emplace_back(vp,N);

        if (nbytes < m_hunk)
        {
            //
            // Add leftover chunk to free list.
            //
            // Insert with a hint -- should be largest block in the set.
            //
            void* block = static_cast<char*>(vp) + nbytes;

            m_freelist.insert(m_freelist.end(), Node(block, vp, m_hunk-nbytes));
        }

        m_busylist.insert(Node(vp, vp, nbytes, stat));
    }
    else
    {
        BL_ASSERT((*free_it).size() >= nbytes);
        BL_ASSERT(m_busylist.find(*free_it) == m_busylist.end());

        vp = (*free_it).block();
        m_busylist.insert(Node(vp, free_it->owner(), nbytes, stat));

        if ((*free_it).size() > nbytes)
        {
            //
            // Insert remainder of free block back into freelist.
            //
            // Insert with a hint -- right after the current block being split.
            //
            Node freeblock = *free_it;

            freeblock.size(freeblock.size() - nbytes);

            freeblock.block(static_cast<char*>(vp) + nbytes);

            m_freelist.insert(free_it, freeblock);
        }

        m_freelist.erase(free_it);
    }

    m_actually_used += nbytes;

    BL_ASSERT(vp != nullptr);

    return vp;
}

std::pair<void*,std::size_t>
CArena::alloc_in_place (void* pt, std::size_t szmin, std::size_t szmax)
{
    std::lock_guard<std::mutex> lock(carena_mutex);

    std::size_t nbytes_max = Arena::align(szmax == 0 ? 1 : szmax);

    if (pt != nullptr) { // Try to allocate in-place first
        auto busy_it = m_busylist.find(Node(pt,nullptr,0));
        if (busy_it == m_busylist.end()) {
            amrex::Abort("CArena::alloc_in_place: unknown pointer");
            return std::make_pair(nullptr,0);
        }
        AMREX_ASSERT(m_freelist.find(*busy_it) == m_freelist.end());

        if (busy_it->size() >= szmax) {
            return std::make_pair(pt, busy_it->size());
        }

        void* next_block = (char*)pt + busy_it->size();
        auto next_it = m_freelist.find(Node(next_block,nullptr,0));
        if (next_it != m_freelist.end() && busy_it->coalescable(*next_it)) {
            std::size_t total_size = busy_it->size() + next_it->size();
            if (total_size >= szmax) {
                // Must use nbytes_max instead of szmax for alignment.
                std::size_t new_size = std::min(total_size, nbytes_max);
                std::size_t left_size = total_size - new_size;
                if (left_size <= 64) {
                    m_freelist.erase(next_it);
                    new_size = total_size;
                } else {
                    auto& free_node = const_cast<Node&>(*next_it);
                    free_node.block((char*)pt + new_size);
                    free_node.size(left_size);
                }
#ifdef AMREX_TINY_PROFILING
                if (m_do_profiling) {
                    TinyProfiler::memory_free(busy_it->size(), busy_it->mem_stat());
                    auto* stat = TinyProfiler::memory_alloc(new_size,
                                                            m_profiling_stats);
                    const_cast<Node&>(*busy_it).mem_stat(stat);
                }
#endif
                m_actually_used += new_size - busy_it->size();
                const_cast<Node&>(*busy_it).size(new_size);
                return std::make_pair(pt, new_size);
            } else if (total_size >= szmin) {
                m_freelist.erase(next_it);
#ifdef AMREX_TINY_PROFILING
                if (m_do_profiling) {
                    TinyProfiler::memory_free(busy_it->size(), busy_it->mem_stat());
                    auto* stat = TinyProfiler::memory_alloc(total_size,
                                                            m_profiling_stats);
                    const_cast<Node&>(*busy_it).mem_stat(stat);
                }
#endif
                m_actually_used += total_size - busy_it->size();
                const_cast<Node&>(*busy_it).size(total_size);
                return std::make_pair(pt, total_size);
            }
        }

        if (busy_it->size() >= szmin) {
            return std::make_pair(pt, busy_it->size());
        }
    }

    void* newp = alloc_protected(nbytes_max);
    return std::make_pair(newp, nbytes_max);
}

void*
CArena::shrink_in_place (void* pt, std::size_t new_size)
{
    if ((pt == nullptr) || (new_size == 0)) { return nullptr; }

    new_size = Arena::align(new_size);

    std::lock_guard<std::mutex> lock(carena_mutex);

    auto busy_it = m_busylist.find(Node(pt,nullptr,0));
    if (busy_it == m_busylist.end()) {
        amrex::Abort("CArena::shrink_in_place: unknown pointer");
        return nullptr;
    }
    AMREX_ASSERT(m_freelist.find(*busy_it) == m_freelist.end());

    auto const old_size = busy_it->size();

    if (new_size > old_size) {
        amrex::Abort("CArena::shrink_in_place: wrong size. Cannot shrink to a larger size.");
        return nullptr;
    } else if (new_size == old_size) {
        return pt;
    } else {
        auto const leftover_size = old_size - new_size;

        void* pt2 = static_cast<char*>(pt) + new_size;
        Node new_free_node(pt2, busy_it->owner(), leftover_size);

        void* pt_end = static_cast<char*>(pt) + old_size;
        auto free_it = m_freelist.find(Node(pt_end,nullptr,0));
        if ((free_it == m_freelist.end()) || ! new_free_node.coalescable(*free_it)) {
            m_freelist.insert(free_it, new_free_node);
        } else {
            auto& node = const_cast<Node&>(*free_it);
            // This is safe because the free list is std::set and the
            // modification of `block` does not change the order of elements
            // in the container, even though Node's operator< uses block.
            node.block(pt2);
            node.size(leftover_size + node.size());
        }

        const_cast<Node&>(*busy_it).size(new_size);

        m_actually_used -= leftover_size;

#ifdef AMREX_TINY_PROFILING
        if (m_do_profiling) {
            TinyProfiler::memory_free(old_size, busy_it->mem_stat());
            auto* stat = TinyProfiler::memory_alloc(new_size, m_profiling_stats);
            const_cast<Node&>(*busy_it).mem_stat(stat);
        }
#endif

        return pt;
    }
}

void
CArena::free (void* vp)
{
    if (vp == nullptr) {
        //
        // Allow calls with NULL as allowed by C++ delete.
        //
        return;
    }

    std::lock_guard<std::mutex> lock(carena_mutex);

    //
    // `vp' had better be in the busy list.
    //
    auto busy_it = m_busylist.find(Node(vp,nullptr,0));
    if (busy_it == m_busylist.end()) {
        amrex::Abort("CArena::free: unknown pointer");
        return;
    }
    BL_ASSERT(m_freelist.find(*busy_it) == m_freelist.end());

    m_actually_used -= busy_it->size();

#ifdef AMREX_TINY_PROFILING
    TinyProfiler::memory_free(busy_it->size(), busy_it->mem_stat());
#endif

    //
    // Put free'd block on free list and save iterator to insert()ed position.
    //
    std::pair<NL::iterator,bool> pair_it = m_freelist.insert(*busy_it);

    BL_ASSERT(pair_it.second == true);

    auto free_it = pair_it.first;

    BL_ASSERT(free_it != m_freelist.end() && (*free_it).block() == (*busy_it).block());
    //
    // And remove from busy list.
    //
    m_busylist.erase(busy_it);
    //
    // Coalesce freeblock(s) on lo and hi side of this block.
    //
    if (!(free_it == m_freelist.begin()))
    {
        auto lo_it = free_it;

        --lo_it;

        void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

        if (addr == (*free_it).block() && lo_it->coalescable(*free_it))
        {
            //
            // This cast is needed as iterators to set return const values;
            // i.e. we can't legally change an element of a set.
            // In this case I want to change the size() of a block
            // in the freelist.  Since size() is not used in the ordering
            // relations in the set, this won't effect the order;
            // i.e. it won't muck up the ordering of elements in the set.
            // I don't want to have to remove the element from the set and
            // then reinsert it with a different size() as it'll just go
            // back into the same place in the set.
            //
            Node* node = const_cast<Node*>(&(*lo_it));
            BL_ASSERT(node != nullptr);
            node->size((*lo_it).size() + (*free_it).size());
            m_freelist.erase(free_it);
            free_it = lo_it;
        }
    }

    auto hi_it = free_it;

    void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

    if (++hi_it != m_freelist.end() && addr == (*hi_it).block() && hi_it->coalescable(*free_it))
    {
        //
        // Ditto the above comment.
        //
        Node* node = const_cast<Node*>(&(*free_it));
        BL_ASSERT(node != nullptr);
        node->size((*free_it).size() + (*hi_it).size());
        m_freelist.erase(hi_it);
    }
}

std::size_t
CArena::freeUnused ()
{
    std::lock_guard<std::mutex> lock(carena_mutex);
    return freeUnused_protected();
}

std::size_t
CArena::freeUnused_protected ()
{
    std::size_t nbytes = 0;
    m_alloc.erase(std::remove_if(m_alloc.begin(), m_alloc.end(),
                                 [&nbytes,this] (std::pair<void*,std::size_t> a)
                                 {
                                     // We cannot simply use std::set::erase because
                                     // Node::operator== only compares the starting address.
                                     auto it = m_freelist.find(Node(a.first,nullptr,0));
                                     if (it != m_freelist.end() &&
                                         it->owner() == a.first &&
                                         it->size()  == a.second)
                                     {
                                         it = m_freelist.erase(it);
                                         nbytes += a.second;
                                         deallocate_system(a.first,a.second);
                                         return true;
                                     }
                                     return false;
                                 }),
                  m_alloc.end());
    m_used -= nbytes;
    return nbytes;
}

bool
CArena::hasFreeDeviceMemory (std::size_t sz)
{
#ifdef AMREX_USE_GPU
    if (isDevice() || isManaged()) {
        std::lock_guard<std::mutex> lock(carena_mutex);

        std::size_t nbytes = Arena::align(sz == 0 ? 1 : sz);

        if (static_cast<Long>(m_used+nbytes) >= arena_info.release_threshold) {
            freeUnused_protected();
        }

        //
        // Find node in freelist at lowest memory address that'll satisfy request.
        //
        NL::iterator free_it = m_freelist.begin();

        for ( ; free_it != m_freelist.end(); ++free_it) {
            if ((*free_it).size() >= nbytes) {
                break;
            }
        }

        if (free_it == m_freelist.end()) {
            const std::size_t N = nbytes < m_hunk ? m_hunk : nbytes;
            return Gpu::Device::freeMemAvailable() > N;
        } else {
            return true;
        }
    } else
#endif
    {
        amrex::ignore_unused(sz);
        return true;
    }
}

void
CArena::registerForProfiling ([[maybe_unused]] const std::string& memory_name)
{
#ifdef AMREX_TINY_PROFILING
    m_do_profiling = true;
    TinyProfiler::RegisterArena(memory_name, m_profiling_stats);
#endif
}

std::size_t
CArena::heap_space_used () const noexcept
{
    return m_used;
}

std::size_t
CArena::heap_space_actually_used () const noexcept
{
    return m_actually_used;
}

std::size_t
CArena::sizeOf (void* p) const noexcept
{
    if (p == nullptr) {
        return 0;
    } else {
        auto it = m_busylist.find(Node(p,nullptr,0));
        if (it == m_busylist.end()) {
            return 0;
        } else {
            return it->size();
        }
    }
}

void
CArena::PrintUsage (std::string const& name) const
{
    Long min_megabytes = static_cast<Long>(heap_space_used() / (1024*1024));
    Long max_megabytes = min_megabytes;
    Long actual_min_megabytes = static_cast<Long>(heap_space_actually_used() / (1024*1024));
    Long actual_max_megabytes = actual_min_megabytes;
    const int IOProc = ParallelDescriptor::IOProcessorNumber();
    ParallelReduce::Min<Long>({min_megabytes, actual_min_megabytes},
                              IOProc, ParallelDescriptor::Communicator());
    ParallelReduce::Max<Long>({max_megabytes, actual_max_megabytes},
                              IOProc, ParallelDescriptor::Communicator());
#ifdef AMREX_USE_MPI
    amrex::Print() << "[" << name << "] space (MB) allocated spread across MPI: ["
                   << min_megabytes << " ... " << max_megabytes << "]\n"
                   << "[" << name << "] space (MB) used      spread across MPI: ["
                   << actual_min_megabytes << " ... " << actual_max_megabytes << "]\n";
#else
    amrex::Print() << "[" << name << "] space allocated (MB): " << min_megabytes << "\n";
    amrex::Print() << "[" << name << "] space used      (MB): " << actual_min_megabytes << "\n";
#endif
}

void
CArena::PrintUsage (std::ostream& os, std::string const& name, std::string const& space) const
{
    auto megabytes = heap_space_used() / (1024*1024);
    auto actual_megabytes = heap_space_actually_used() / (1024*1024);
    os << space << "[" << name << "] space allocated (MB): " << megabytes << "\n";
    os << space << "[" << name << "] space used      (MB): " << actual_megabytes << "\n";
    os << space << "[" << name << "]: " << m_alloc.size() << " allocs, "
       << m_busylist.size() << " busy blocks, " << m_freelist.size() << " free blocks\n";
}

std::ostream& operator<< (std::ostream& os, const CArena& arena)
{
    os << "CArea:\n"
       << "    Hunk size: " << arena.m_hunk << "\n"
       << "    Memory allocated: " << arena.m_used << "\n"
       << "    Memory actually used: " << arena.m_actually_used << "\n";

    if (arena.m_alloc.empty()) {
        os << "    No memory allocations\n";
    } else {
        os << "    List of memory alloations: (address, size)\n";
        for (auto const& a : arena.m_alloc) {
            os << "        " << a.first << ", " << a.second << "\n";
        }
    }

    if (arena.m_freelist.empty()) {
        os << "    No free nodes\n";
    } else {
        os << "    List of free nodes: (address, owner, size)\n";
        for (auto const& a : arena.m_freelist) {
            os << "        " << a.block() << ", " << a.owner() << ", "
               << a.size() << "\n";
        }
    }

    if (arena.m_busylist.empty()) {
        os << "    No busy nodes\n";
    } else {
        os << "    List of busy nodes: (address, owner, size)\n";
        for (auto const& a : arena.m_busylist) {
            os << "        " << a.block() << ", " << a.owner() << ", "
               << a.size() << "\n";
        }
    }

    return os;
}

}
