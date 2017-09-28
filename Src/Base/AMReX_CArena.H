#ifndef BL_CARENA_H
#define BL_CARENA_H

#include <cstddef>
#include <set>
#include <vector>

#include <AMReX_Arena.H>
#include <AMReX_Array.H>

namespace amrex {

/**
* \brief A Concrete Class for Dynamic Memory Management
* This is a coalescing memory manager.  It allocates (possibly) large
* chunks of heap space and apportions it out as requested.  It merges
* together neighboring chunks on each free().
*/

class CArena
    :
    public Arena
{
public:
    /**
    * \brief Construct a coalescing memory manager.  hunk_size is the
    * minimum size of hunks of memory to allocate from the heap.
    * If hunk_size == 0 we use DefaultHunkSize as specified below.
    */
    CArena (size_t hunk_size = 0);
    
    //! The destructor.
    virtual ~CArena () override;

    //! Allocate some memory.
    virtual void* alloc (size_t nbytes) override;

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free (void* ap) override;

    //! The current amount of heap space used by the CArena object.
    size_t heap_space_used () const;

    //! The default memory hunk size to grab from the heap.
    enum { DefaultHunkSize = 1024*1024*8 };

protected:
    //! The nodes in our free list and block list.
    class Node
    {
    public:
        //! The default constructor.
        Node ()
            :
            m_block(0), m_size(0) {}
        
        //! Another constructor.
        Node (void* a_block, size_t a_size)
            :
            m_block(a_block), m_size(a_size) {}

        //! The copy constructor.
        Node (const Node& rhs)
            :
            m_block(rhs.m_block), m_size(rhs.m_size) {}

        //! The copy assignment constructor.
        Node& operator= (const Node& rhs)
        {
            m_block = rhs.m_block;
            m_size  = rhs.m_size;
            return *this;
        }

        //! The "less-than" operator.
        bool operator< (const Node& rhs) const
        {
            return m_block < rhs.m_block;
        }

        //! The equality operator. 
        bool operator== (const Node& rhs) const
        {
            return m_block == rhs.m_block;
        }

        //! The block address.
        void* block () const { return m_block; }

        //! Set block address.
        void block (void* blk) { m_block = blk; }

        //! The size of the memory block.
        size_t size () const { return m_size; }

        //! Set size.
        void size (size_t sz) { m_size = sz; }

    private:
        //! The block of memory we reference.
        void* m_block;
        //! The size of the block we represent.
        size_t m_size;
    };

    /**
    * \brief The type of our freelist and blocklist.
    * We use a set sorted from lo to hi memory addresses.
    */
    typedef std::set<Node> NL;

    //! The list of blocks allocated via ::operator new().
    std::vector<void*> m_alloc;

    /**
    * \brief The free list of allocated but not currently used blocks.
    * Maintained in lo to hi memory sorted order.
    */
    NL m_freelist;

    /**
    * \brief The list of busy blocks.
    * A block is either on the freelist or on the blocklist, but not on both.
    */
    NL m_busylist;
    //! The minimal size of hunks to request via ::operator new().
    size_t m_hunk;
    //! The amount of heap space currently allocated.
    size_t m_used;

private:
    //! Disallowed.
    CArena (const CArena& rhs);
    CArena& operator= (const CArena& rhs);
};

}

#endif /*BL_CARENA_H*/
