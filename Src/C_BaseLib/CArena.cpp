//BL_COPYRIGHT_NOTICE

//
// $Id: CArena.cpp,v 1.6 1998-02-08 20:03:59 lijewski Exp $
//

#ifdef BL_USE_NEW_HFILES
#include <iostream>
#include <iomanip>
#include <fstream>
using std::ofstream;
using std::ios;
using std::setw;
using std::clog;
using std::endl;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#endif

#include <CArena.H>

static CArena The_Static_FAB_CArena;

extern CArena* The_FAB_CArena = &The_Static_FAB_CArena;

CArena::CArena (size_t hunk_size)
{
    //
    // Force alignment of hunksize.
    //
    m_hunk = Arena::align(hunk_size == 0 ? DefaultHunkSize : hunk_size);

    assert(m_hunk >= hunk_size);
    assert(m_hunk%sizeof(Arena::Word) == 0);
}

CArena::~CArena ()
{
    for (int i = 0; i < m_alloc.size(); i++)
    {
        ::operator delete(m_alloc[i]);
    }
}

void*
CArena::alloc (size_t nbytes,
               void** client)
{
    assert(!(client == 0));

    nbytes = Arena::align(nbytes);

    NL::iterator free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it)
    {
        if ((*free_it).size() >= nbytes)
            break;
    }

    if (free_it == m_freelist.end())
    {
        void* vp = ::operator new(nbytes < m_hunk ? m_hunk : nbytes);

        m_alloc.push_back(vp);

        *client = vp;

        m_busylist.push_back(Node(vp, client, nbytes));

        if (nbytes < m_hunk)
        {
            //
            // Add leftover chunk to free list.
            //
            void* block = static_cast<char*>(vp) + nbytes;

            m_freelist.push_front(Node(block, 0, m_hunk - nbytes));
        }
    }
    else
    {
        assert((*free_it).client() == 0);

        *client = (*free_it).block();

        m_busylist.push_back(Node(*client, client, nbytes));

        if ((*free_it).size() == nbytes)
        {
            m_freelist.erase(free_it);
        }
        else
        {
            //
            // Update freelist block to reflect space unused by busyblock.
            //
            (*free_it).size((*free_it).size() - nbytes);

            (*free_it).block(static_cast<char*>(*client) + nbytes);
        }
    }

    assert(!(*client == 0));

    return *client;
}

void
CArena::free (void* vp)
{
    if (vp == 0)
        //
        // Allow calls with NULL as allowed by C++ delete.
        //
        return;

    NL::iterator busy_it = m_busylist.begin();

    for ( ; busy_it != m_busylist.end(); ++busy_it)
    {
        if ((*busy_it).block() == vp)
            break;
    }
    assert(!(busy_it == m_busylist.end()));

    void* free_block = static_cast<char*>((*busy_it).block());

    m_freelist.push_front(Node(free_block, 0, (*busy_it).size()));

    m_busylist.erase(busy_it);
    //
    // Coalesce freeblock(s) on lo and hi side of this block.
    //
    // Sort from lo to hi memory addresses.
    //
    m_freelist.sort();

    NL::iterator free_it = m_freelist.begin();

    for ( ; free_it != m_freelist.end(); ++free_it)
    {
        if ((*free_it).block() == free_block)
            break;
    }
    assert(!(free_it == m_freelist.end()));

    if (free_it != m_freelist.begin())
    {
        NL::iterator lo_it = free_it;

        --lo_it;

        void* addr = static_cast<char*>((*lo_it).block()) + (*lo_it).size();

        if (addr == (*free_it).block())
        {
            (*lo_it).size((*lo_it).size() + (*free_it).size());

            m_freelist.erase(free_it);

            free_it = lo_it;
        }
    }

    NL::iterator hi_it = free_it;

    void* addr = static_cast<char*>((*free_it).block()) + (*free_it).size();

    if (++hi_it != m_freelist.end() && addr == (*hi_it).block())
    {
        (*free_it).size((*free_it).size() + (*hi_it).size());

        m_freelist.erase(hi_it);
    }
}

void
CArena::compact ()
{

}
