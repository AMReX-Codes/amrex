//BL_COPYRIGHT_NOTICE

#ifndef BL_CARENA
#define BL_CARENA

#ifdef BL_USE_NEW_HFILES
#include <cstddef>
#include <iostream>
#include <list>
#include <vector>
using std::list;
using std::ostream;
using std::vector;
#else
#include <stddef.h>
#include <iostream.h>
#include <list.h>
#include <vector.h>
#endif

#include <Boolean.H>
#include <Arena.H>
#include <Array.H>

class CArena
    :
    public Arena
{
public:
    //
    // `hunk_size' is the minimum size of hunks of memory to
    // allocate in any call of alloc().  If hunk_size == 0 we use
    // DefaultHunkSize as specified below.
    //
    CArena (size_t hunk_size = 0);
    //
    // The destructor.
    //
    ~CArena ();
    //
    // Allocate some memory.
    //
    // `ind' is the address of the memory in the user's code.
    // We use this when we compress memory.
    // If user's cache data pointers, they're (sooner or later) hosed.
    //
    virtual void* alloc (size_t nbytes, void** ind);
    //
    // Free up allocated memory.
    //
    virtual void free (void* ap);
    //
    // Compact the `used' memory.
    //
    void compact ();

protected:
    //
    // If hunk_size==0 on construction, use this size.
    //
    enum { DefaultHunkSize = 1024*1024 };
    //
    // The nodes in our free list and block list.
    //
    class Node
    {
    public:
        //
        // The default constructor.
        //
        Node ()
            :
            m_block(0), m_client(0), m_size(0) {}
        //
        // Another constructor.
        //
        Node (void* block, void** client, size_t size)
            :
            m_block(block), m_client(client), m_size(size) {}
        //
        // The "less-than" operator.
        //
        bool operator< (const Node& rhs) const
        {
            return m_block < rhs.m_block;
        }
        //
        // The equality operator. 
        //
        bool operator== (const Node& rhs) const
        {
            return m_block == rhs.m_block;
        }
        //
        // The block address.
        //
        void* block () const { return m_block; }
        //
        // Set block address.
        //
        void block (void* blk) { m_block = blk; }
        //
        // The client address.
        //
        void** client () const { return m_client; }
        //
        // Set client address.
        //
        void client (void** clnt) { m_client = clnt; }
        //
        // The size of the memory block.
        //
        size_t size () const { return m_size; }
        //
        // Set size.
        //
        void size (size_t sz) { m_size = sz; }

    private:
        //
        // The block of memory we reference.
        //
        void* m_block;
        //
        // The address of the client pointer that points to m_block. 
        //
        void** m_client;
        //
        // The size of the block we represent.
        //
        size_t m_size;
    };
    //
    // The type of our freelist and blocklist.
    //
    typedef list < Node > NL;
    //
    // The list of blocks allocated via ::operator new().
    //
    vector<void*> m_alloc;
    //
    // The free list of allocated but not currently used blocks.
    //
    NL m_freelist;
    //
    // The list of busy blocks.
    // A block is either on the freelist or on the blocklist, but not on both.
    //
    NL m_busylist;
    //
    // The minimal size of hunks to request via ::operator new().
    //
    size_t m_hunk;

private:
    //
    // Disallowed.
    //
    CArena (const CArena& rhs);
    CArena& operator= (const CArena& rhs);
};

//
// The CArena used by BaseFab code.
//
extern CArena* The_FAB_CArena;

#endif /*BL_CARENA*/
