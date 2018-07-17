
#include <utility>
#include <cstring>

#include <AMReX_CArena.H>
#include <AMReX_BLassert.H>

#ifdef AMREX_USE_DEVICE
#include <AMReX_Device.H>
#endif

namespace amrex {

CArena::CArena (size_t hunk_size)
{
    //
    // Force alignment of hunksize.
    //
    m_hunk = Arena::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);
    m_used = 0;

    BL_ASSERT(m_hunk >= hunk_size);
    BL_ASSERT(m_hunk%Arena::align_size == 0);
}

CArena::~CArena ()
{
    for (unsigned int i = 0, N = m_alloc.size(); i < N; i++)
#ifdef AMREX_USE_CUDA
	if (device_use_hostalloc)
	    gpu_freehost(m_alloc[i]);
	else
	    gpu_free(m_alloc[i]);
#else
        ::operator delete(m_alloc[i]);
#endif
}

void*
CArena::alloc (size_t nbytes)
{
    nbytes = Arena::align(nbytes == 0 ? 1 : nbytes);
    //
    // Find node in freelist at lowest memory address that'll satisfy request.
    //
    NL::iterator free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it)
        if ((*free_it).size() >= nbytes)
            break;

    void* vp = 0;

    if (free_it == m_freelist.end())
    {
        const size_t N = nbytes < m_hunk ? m_hunk : nbytes;

#if (defined(AMREX_USE_CUDA) && defined(AMREX_USE_CUDA_UM))
        if (device_use_hostalloc) {

	    gpu_hostalloc(&vp, &N);

	}
	else if (device_use_managed_memory) {

	    gpu_malloc_managed(&vp, &N);
	    if (device_set_readonly)
		Device::mem_advise_set_readonly(vp, N);
	    if (device_set_preferred) {
		const int device = Device::deviceId();
		Device::mem_advise_set_preferred(vp, N, device);
	    }

	}
	else {

	    gpu_malloc(&vp, &N);

	}
#else
        vp = ::operator new(N);
#endif

        m_used += N;

        m_alloc.push_back(vp);

        if (nbytes < m_hunk)
        {
            //
            // Add leftover chunk to free list.
            //
            // Insert with a hint -- should be largest block in the set.
            //
            void* block = static_cast<char*>(vp) + nbytes;

            m_freelist.insert(m_freelist.end(), Node(block, m_hunk-nbytes));
        }
    }
    else
    {
        BL_ASSERT((*free_it).size() >= nbytes);
        BL_ASSERT(m_busylist.find(*free_it) == m_busylist.end());

        vp = (*free_it).block();

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

    m_busylist.insert(Node(vp, nbytes));

    BL_ASSERT(!(vp == 0));

    return vp;
}

void
CArena::free (void* vp)
{
    if (vp == 0)
        //
        // Allow calls with NULL as allowed by C++ delete.
        //
        return;
    //
    // `vp' had better be in the busy list.
    //
    NL::iterator busy_it = m_busylist.find(Node(vp,0));

    BL_ASSERT(!(busy_it == m_busylist.end()));
    BL_ASSERT(m_freelist.find(*busy_it) == m_freelist.end());
    //
    // Put free'd block on free list and save iterator to insert()ed position.
    //
    std::pair<NL::iterator,bool> pair_it = m_freelist.insert(*busy_it);

    BL_ASSERT(pair_it.second == true);

    NL::iterator free_it = pair_it.first;

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
        NL::iterator lo_it = free_it;

        --lo_it;

        void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

        if (addr == (*free_it).block())
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
            BL_ASSERT(!(node == 0));
            node->size((*lo_it).size() + (*free_it).size());
            m_freelist.erase(free_it);
            free_it = lo_it;
        }
    }

    NL::iterator hi_it = free_it;

    void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

    if (++hi_it != m_freelist.end() && addr == (*hi_it).block())
    {
        //
        // Ditto the above comment.
        //
        Node* node = const_cast<Node*>(&(*free_it));
        BL_ASSERT(!(node == 0));
        node->size((*free_it).size() + (*hi_it).size());
        m_freelist.erase(hi_it);
    }
}

#ifdef AMREX_USE_DEVICE
// Device allocators are not currently implemented in CArena.

void*
CArena::alloc_device (size_t nbytes)
{
    void* pt = 0;
    return pt;
}

void
CArena::free_device (void* pt)
{
}
#endif

size_t
CArena::heap_space_used () const
{
    return m_used;
}

}
