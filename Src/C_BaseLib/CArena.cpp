//BL_COPYRIGHT_NOTICE

//
// $Id: CArena.cpp,v 1.5 1998-02-08 18:19:28 lijewski Exp $
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

    NL::iterator it = m_freelist.begin();

    for ( ; it != m_freelist.end(); ++it)
    {
        if ((*it).size() >= nbytes)
            break;
    }

    if (it == m_freelist.end())
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
        assert((*it).client() == 0);

        *client = (*it).block();

        m_busylist.push_back(Node(*client, client, nbytes));

        if ((*it).size() == nbytes)
        {
            m_freelist.erase(it);
        }
        else
        {
            //
            // Update freelist block to reflect space unused by busyblock.
            //
            (*it).size((*it).size() - nbytes);

            (*it).block(static_cast<char*>(*client) + nbytes);
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
    //
    // Block had better be in busylist.
    //
    NL::iterator it = m_busylist.begin();

    for ( ; it != m_busylist.end(); ++it)
    {
        if ((*it).block() == vp)
            break;
    }

    assert(!(it == m_busylist.end()));

    m_freelist.push_front(Node((*it).block(), 0, (*it).size()));
    //
    // Set client pointer to NULL -- try to catch client pointer caching :-)
    //
    //*(*it).client() = 0;

    m_busylist.erase(it);
}

void
CArena::compact ()
{

}
