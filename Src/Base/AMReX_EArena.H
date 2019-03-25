#ifndef AMREX_EARENA_H_
#define AMREX_EARENA_H_

#include <cstddef>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>
#include <mutex>
#include <unordered_set>
#include <functional>

#include <AMReX_Arena.H>

namespace amrex {

/**
* \brief A Concrete Class for Dynamic Memory Management using best fit. 
* This is a coalescing memory manager.  It allocates (possibly) large
* chunks of heap space and apportions it out as requested.  It merges
* together neighboring chunks on each free().
*/

class EArena
    :
    public Arena
{
public:

    EArena (std::size_t hunk_size = 0, ArenaInfo info = ArenaInfo());
    EArena (const EArena& rhs) = delete;
    EArena& operator= (const EArena& rhs) = delete;
    virtual ~EArena () override;

    virtual void* alloc (std::size_t nbytes) override final;
    virtual void free (void* vp) override final;

    //! The current amount of heap space used by the EArena object.
    std::size_t heap_space_used () const noexcept;
    //! Free space available in the arena
    std::size_t free_space_available () const noexcept;

    //! The default memory hunk size to grab from the heap.
    enum { DefaultHunkSize = 1024*1024*8 };

protected:

    struct Node
    {
        Node (void* block, void* owner, std::size_t size)
            : m_block(reinterpret_cast<uintptr_t>(block)),
              m_owner(reinterpret_cast<uintptr_t>(owner)),
              m_size(size)
            {}

        Node (std::size_t size)
            : m_block(0), m_owner(0), m_size(size)
            {}

        uintptr_t m_block;
        uintptr_t m_owner;
        std::size_t m_size;

        struct CompareSize {
            bool operator () (Node const& lhs, Node const& rhs) const noexcept {
                return (lhs.m_size < rhs.m_size)
                    || ((lhs.m_size == rhs.m_size)
                        && (lhs.m_block < rhs.m_block));
            }
        };

        struct CompareAddr {
            bool operator () (Node const& lhs, Node const& rhs) const noexcept {
                return (lhs.m_owner < rhs.m_owner)
                    || ((lhs.m_owner == rhs.m_owner)
                        && (lhs.m_block < rhs.m_block));
            }
        };

        struct hash {
            std::size_t operator() (Node const& n) const noexcept {
                return std::hash<uintptr_t>{}(n.m_block);
            }
        };

        struct equal {
            bool operator() (Node const& lhs, Node const& rhs) const noexcept {
                return lhs.m_block == rhs.m_block;
            }
        };
    };

    std::vector<void*> m_alloc;

    using FreeList = std::set<Node,Node::CompareSize>;
    FreeList m_freelist;

    // This store the same nodes as FreeList, but in different order
    using MergeList = std::set<Node,Node::CompareAddr>;
    MergeList m_mergelist;

    using BusyList = std::unordered_set<Node,Node::hash,Node::equal>;
    BusyList m_busylist;

    //! The minimal size of hunks to request from system
    std::size_t m_hunk;
    //! The amount of heap space currently allocated from system.
    std::size_t m_used_size;
    //! The amount of free space in arena
    std::size_t m_free_size;

    std::mutex earena_mutex;
};

}

#endif
