#ifndef BL_CARENA_H
#define BL_CARENA_H

#include <cstddef>
#include <set>
#include <vector>
#include <mutex>
#include <unordered_set>
#include <functional>

#include <AMReX_Arena.H>

namespace amrex {

/**
* \brief A Concrete Class for Dynamic Memory Management using first fit.
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
    CArena (std::size_t hunk_size = 0, ArenaInfo info = ArenaInfo());

    CArena (const CArena& rhs) = delete;
    CArena& operator= (const CArena& rhs) = delete;
    
    //! The destructor.
    virtual ~CArena () override;

    //! Allocate some memory.
    virtual void* alloc (std::size_t nbytes) override final;

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    virtual void free (void* ap) override final;

    //! The current amount of heap space used by the CArena object.
    std::size_t heap_space_used () const noexcept;

    //! The default memory hunk size to grab from the heap.
    enum { DefaultHunkSize = 1024*1024*8 };

protected:
    //! The nodes in our free list and block list.
    class Node
    {
    public:
        Node (void* a_block, void* a_owner, std::size_t a_size) noexcept
            :
            m_block(a_block), m_owner(a_owner), m_size(a_size) {}

        //! The "less-than" operator.
        bool operator< (const Node& rhs) const noexcept
        {
            return m_block < rhs.m_block;
        }

        //! The equality operator. 
        bool operator== (const Node& rhs) const noexcept
        {
            return m_block == rhs.m_block;
        }

        //! The block address.
        void* block () const noexcept { return m_block; }

        //! Set block address.
        void block (void* blk) noexcept { m_block = blk; }

        //! The size of the memory block.
        std::size_t size () const noexcept { return m_size; }

        //! Set size.
        void size (std::size_t sz) noexcept { m_size = sz; }

        void* owner () const noexcept { return m_owner; }

        bool coalescable (const Node& rhs) const noexcept {
            return m_owner == rhs.m_owner;
        }

        struct hash {
            std::size_t operator() (const Node& n) const noexcept {
                return std::hash<void*>{}(n.m_block);
            }
        };

    private:
        //! The block of memory we reference.
        void* m_block;
        //! The starting address of the original allocation
        void* m_owner;
        //! The size of the block we represent.
        std::size_t m_size;
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
//    NL m_busylist;
    std::unordered_set<Node, Node::hash> m_busylist;
    //! The minimal size of hunks to request from system
    std::size_t m_hunk;
    //! The amount of heap space currently allocated.
    std::size_t m_used;

    std::mutex carena_mutex;
};

}

#endif /*BL_CARENA_H*/
