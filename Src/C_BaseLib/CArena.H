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

#include <Arena.H>
#include <Array.H>

//
//@Man:
//@Memo: A Concrete Class for Dynamic Memory Management
/*@Doc:

  This is a coalescing memory manager.  It allocates (possibly) large
  chunks of heap space and apportions it out as requested.  It merges
  together neighboring chunks on each free().
*/

class CArena
    :
    public Arena
{
public:

    /*@ManDoc: Construct a coalescing memory manager.  `hunk\_size' is the
               minimum size of hunks of memory to allocate from the heap.
               If hunk\_size == 0 we use DefaultHunkSize as specified below.
    */
    CArena (size_t hunk_size = 0);
    //
    //@ManDoc: The destructor.
    //
    ~CArena ();
    //
    //@ManDoc: Allocate some memory.
    //
    virtual void* alloc (size_t nbytes);

    /*@ManDoc: Free up allocated memory.  Merge neighboring free memory chunks
               into largest possible chunk.
    */
    virtual void free (void* ap);

    /*@ManDoc: Mirror the C calloc() function.  Returns zero'd memory of
               size nmemb*size.  This is here so that we can implement
               malloc(3) and family.  User's shouldn't use this function.
    */
    void* calloc (size_t nmemb, size_t size);

    /*@ManDoc: Mirror the C realloc() function.  This is here so that we can
               implement malloc(3) and family.  User's shouldn't use this
               function.
    */
    void* realloc (void* ptr, size_t size);
    //
    //@ManDoc: The default memory hunk size to grab from the heap.
    //
    enum { DefaultHunkSize = 1024*1024 };

protected:
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
