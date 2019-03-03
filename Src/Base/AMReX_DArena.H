#ifndef AMREX_D_ARENA_H_
#define AMREX_D_ARENA_H_

#include <cstddef>
#include <unordered_set>
#include <unordered_map>
#include <array>
#include <mutex>

#include <AMReX_Arena.H>

namespace amrex {

/*
 * \brief Buddy memory allocator
 */

class DArena
    : public Arena
{
public:

    DArena (std::size_t max_size, std::size_t max_block_size, ArenaInfo info);

    DArena (DArena const&) = delete;
    DArena (DArena &&) = delete;
    DArena& operator= (DArena const&) = delete;
    DArena& operator= (DArena&&) = delete;

    virtual ~DArena () override;

    virtual void* alloc (std::size_t nbytes) override final;
    virtual void free (void* p) override final;

    std::size_t totalMem () const { return m_max_size; }
    std::size_t freeMem () const;

private:
    static constexpr int m_max_max_order = 30;
    // buckets of free blocks
    std::array<std::unordered_set<std::ptrdiff_t>, m_max_max_order+1> m_free;
    // map: allocated pointer offset relative to baseptr -> order
    std::unordered_map<std::ptrdiff_t,int> m_used;
    char* m_baseptr;
    std::size_t m_max_size;
    std::size_t m_block_size;
    int m_max_order;
    std::mutex m_mutex;
    bool warning_printed = false;

    std::ptrdiff_t allocate_order (int order);
    void deallocate_order (int order, std::ptrdiff_t offset);
};

}

#endif
