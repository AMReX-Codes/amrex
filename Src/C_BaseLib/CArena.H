//BL_COPYRIGHT_NOTICE

#ifndef BL_CARENA
#define BL_CARENA

#ifdef BL_USE_NEW_HFILES
#ifdef WIN32
#pragma warning (disable : 4786)
#endif
#include <cstddef>
#include <iostream>
#include <set>
#include <vector>
using std::set;
using std::ostream;
using std::vector;
#else
#include <stddef.h>
#include <iostream.h>
#include <set.h>
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
    virtual void* alloc (size_t nbytes);
    //
    // Free up allocated memory.
    //
    virtual void free (void* ap);

protected:

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
            m_block(0), m_size(0) {}
        //
        // Another constructor.
        //
        Node (void* block, size_t size)
            :
            m_block(block), m_size(size) {}
        //
        // The copy constructor.
        //
        Node (const Node& rhs)
            :
            m_block(rhs.m_block), m_size(rhs.m_size) {}
        //
        // The copy assignment constructor.
        //
        Node& operator= (const Node& rhs)
        {
            m_block = rhs.m_block;
            m_size  = rhs.m_size;
            return *this;
        }
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
        // The size of the block we represent.
        //
        size_t m_size;
    };
    //
    // The type of our freelist and blocklist.
    // We use a set sorted from lo to hi memory addresses.
    //
    typedef set < Node > NL;
    //
    // The list of blocks allocated via ::operator new().
    //
    vector<void*> m_alloc;
    //
    // The free list of allocated but not currently used blocks.
    // Maintained in lo to hi memory sorted order.
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
// The Arena used by BaseFab code.
//
extern Arena* The_FAB_Arena;

#endif /*BL_CARENA*/
