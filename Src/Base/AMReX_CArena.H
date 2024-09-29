#ifndef BL_CARENA_H
#define BL_CARENA_H
#include <AMReX_Config.H>

#include <AMReX_Arena.H>

#include <cstddef>
#include <functional>
#include <iosfwd>
#include <map>
#include <mutex>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

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
    CArena (CArena&& rhs) = delete;
    CArena& operator= (const CArena& rhs) = delete;
    CArena& operator= (CArena&& rhs) = delete;

    //! The destructor.
    ~CArena () override;

    //! Allocate some memory.
    [[nodiscard]] void* alloc (std::size_t nbytes) final;

    /**
     * Try to allocate in-place by extending the capacity of given pointer.
     */
    [[nodiscard]] std::pair<void*,std::size_t>
    alloc_in_place (void* pt, std::size_t szmin, std::size_t szmax) final;

    /**
     * Try to shrink in-place
     */
    [[nodiscard]] void*
    shrink_in_place (void* pt, std::size_t new_size) final;

    /**
    * \brief Free up allocated memory.  Merge neighboring free memory chunks
    * into largest possible chunk.
    */
    void free (void* vp) final;

    std::size_t freeUnused () final;

    /**
     * \brief Does the device have enough free memory for allocating this
     * much memory?  For CPU builds, this always return true.  This is not a
     * const function because it may attempt to release memory back to the
     * system.
     */
    [[nodiscard]] bool hasFreeDeviceMemory (std::size_t sz) final;

    //! The current amount of heap space used by the CArena object.
    std::size_t heap_space_used () const noexcept;

    //! Return the total amount of memory given out via alloc.
    std::size_t heap_space_actually_used () const noexcept;

    //! Return the amount of memory in this pointer.  Return 0 for unknown pointer.
    std::size_t sizeOf (void* p) const noexcept;

    void PrintUsage (std::string const& name) const;

    void PrintUsage (std::ostream& os, std::string const& name, std::string const& space) const;

    //! The default memory hunk size to grab from the heap.
    constexpr static std::size_t DefaultHunkSize = 1024*1024*8;

protected:

    void* alloc_protected (std::size_t nbytes);

    std::size_t freeUnused_protected () final;

    //! The nodes in our free list and block list.
    class Node
    {
    public:
        Node (void* a_block, void* a_owner, std::size_t a_size, MemStat* a_stat=nullptr) noexcept
            :
            m_block(a_block), m_owner(a_owner), m_size(a_size), m_stat(a_stat) {}

        //! The "less-than" operator.
        bool operator< (const Node& rhs) const noexcept
        {
            return std::less<>{}(m_block, rhs.m_block);
        }

        //! The equality operator.
        bool operator== (const Node& rhs) const noexcept
        {
            return m_block == rhs.m_block;
        }

        //! The block address.
        [[nodiscard]] void* block () const noexcept { return m_block; }

        //! Set block address.
        void block (void* blk) noexcept { m_block = blk; }

        //! The size of the memory block.
        [[nodiscard]] std::size_t size () const noexcept { return m_size; }

        //! Set size.
        void size (std::size_t sz) noexcept { m_size = sz; }

        [[nodiscard]] void* owner () const noexcept { return m_owner; }

        [[nodiscard]] bool coalescable (const Node& rhs) const noexcept {
            return m_owner == rhs.m_owner;
        }

        //! Get the MemStat object of the function where this block was allocated.
        [[nodiscard]] MemStat* mem_stat () const { return m_stat; }

        //! Set MemStat
        void mem_stat (MemStat* a_stat) noexcept { m_stat = a_stat; }

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
        //! Used for profiling if this Node represents a user allocated block of memory.
        MemStat* m_stat;
    };

    //! The list of blocks allocated via ::operator new().
    std::vector<std::pair<void*,std::size_t> > m_alloc;

    /**
    * \brief The type of our freelist and blocklist.
    * We use a set sorted from lo to hi memory addresses.
    */
    using NL = std::set<Node>;

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
    std::size_t m_used{0};
    //! The amount of memory given out via alloc().
    std::size_t m_actually_used{0};


    std::mutex carena_mutex;

    friend std::ostream& operator<< (std::ostream& os, const CArena& arena);
};

}

#endif /*BL_CARENA_H*/
