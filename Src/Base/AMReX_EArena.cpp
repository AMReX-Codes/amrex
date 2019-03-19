
#include <utility>
#include <cstring>

#include <AMReX_EArena.H>
#include <AMReX_BLassert.H>
#include <AMReX_Gpu.H>

namespace amrex {

EArena::EArena (std::size_t hunk_size, ArenaInfo info)
{
    arena_info = info;
    m_hunk = Arena::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);
    m_used_size = 0;
    m_free_size = 0;

    AMREX_ASSERT(m_hunk >= hunk_size);
    AMREX_ASSERT(m_hunk%Arena::align_size == 0);
}

EArena::~EArena ()
{
    for (unsigned int i = 0, N = m_alloc.size(); i < N; i++) {
        deallocate_system(m_alloc[i]);
    }
}

void*
EArena::alloc (std::size_t nbytes)
{
    std::lock_guard<std::mutex> lock(earena_mutex);

    nbytes = Arena::align(nbytes == 0 ? 1 : nbytes);

    auto fit = m_freelist.lower_bound(Node{nbytes});

    void* vp = nullptr;

    if (fit == m_freelist.end())
    { // We have to allocate from the system
        const std::size_t N = nbytes < m_hunk ? m_hunk : nbytes;
        vp = allocate_system(N);
        m_used_size += N;
        m_alloc.push_back(vp);

        if (nbytes < N) { // add leftover to free list
            void* block = static_cast<char*>(vp) + nbytes;
            m_freelist.emplace(block, vp, N-nbytes);
            m_mergelist.emplace(block, vp, N-nbytes);
            m_free_size += N-nbytes;
        }

        m_busylist.emplace(vp, vp, nbytes);
    }
    else
    { // found free block
        AMREX_ASSERT(fit->m_size >= nbytes);
        AMREX_ASSERT(m_busylist.find(*fit) == m_busylist.end());

        vp = reinterpret_cast<void*>(fit->m_block);
        void* op = reinterpret_cast<void*>(fit->m_owner);
        m_busylist.emplace(vp, op, nbytes);

        if (fit->m_size > nbytes) {
            void* bp = static_cast<char*>(vp) + nbytes;
            m_freelist.emplace(bp,op,fit->m_size-nbytes);
            m_mergelist.emplace(bp,op,fit->m_size-nbytes);
            m_free_size += fit->m_size-nbytes;
        }

        m_free_size -= fit->m_size;
        m_mergelist.erase(*fit);
        m_freelist.erase(fit);
    }

    AMREX_ASSERT(vp != nullptr);
    return vp;
}

void
EArena::free (void* vp)
{
    std::lock_guard<std::mutex> lock(earena_mutex);
    if (vp == nullptr) return;

    auto bit = m_busylist.find(Node{vp,nullptr,0});
    AMREX_ASSERT(bit != m_busylist.end()); // assert pointer is in busy list
    AMREX_ASSERT(m_freelist.find(*bit) == m_freelist.end()); // assert not in free list

    auto pair_fitb = m_freelist.insert(*bit);
    auto pair_mitb = m_mergelist.insert(*bit);
    m_free_size += bit->m_size;
    AMREX_ASSERT(pair_fitb.second == true && pair_mitb.second);
    auto fit = pair_fitb.first;
    auto mit = pair_mitb.first;

    m_busylist.erase(bit);

    // try to merge
    if (mit != m_mergelist.begin())
    {
        auto lo_it = mit;
        --lo_it;
        if (lo_it->m_owner == mit->m_owner)
        {
            void* plo = reinterpret_cast<void*>(lo_it->m_block);
            void* phi = reinterpret_cast<void*>(mit->m_block);
            void* addr = static_cast<char*>(plo) + lo_it->m_size;
            if (addr == phi)
            {
                m_freelist.erase(*lo_it);
                m_freelist.erase(*mit);
                Node* node = const_cast<Node*>(&(*lo_it));
                node->m_size = lo_it->m_size + mit->m_size;
                m_mergelist.erase(mit);
                mit = lo_it;
                m_freelist.insert(*mit);
            }
        }
    }
    
    auto hi_it = mit;
    ++hi_it;
    if (hi_it != m_mergelist.end() && mit->m_owner == hi_it->m_owner)
    {
       void* plo = reinterpret_cast<void*>(mit->m_block);
       void* phi = reinterpret_cast<void*>(hi_it->m_block);
       void* addr = static_cast<char*>(plo) + mit->m_size;
       if (addr == phi)
       {
           m_freelist.erase(*mit);
           m_freelist.erase(*hi_it);
           Node* node = const_cast<Node*>(&(*mit));
           node->m_size = mit->m_size + hi_it->m_size;
           m_mergelist.erase(hi_it);
           m_freelist.insert(*mit);
       }
    }
}

std::size_t
EArena::heap_space_used () const noexcept
{
    return m_used_size;
}

std::size_t
EArena::free_space_available () const noexcept
{
    return m_free_size;
}

}
