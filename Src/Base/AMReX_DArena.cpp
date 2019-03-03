
#include <cstdlib>
#include <cmath>

#include <AMReX_DArena.H>
#include <AMReX_BLassert.H>
#include <AMReX_Gpu.H>
#include <AMReX_Print.H>
#include <AMReX.H>

namespace amrex {

namespace {
    std::ptrdiff_t make_buddy (int order, std::ptrdiff_t offset)
    {
        std::size_t tmp = offset;
        tmp ^= 1ul << order;
        return static_cast<std::ptrdiff_t>(tmp);
    }
}

DArena::DArena (std::size_t max_size, std::size_t max_block_size, ArenaInfo info)
{
    arena_info = info;

    m_max_order = 0;
    while (max_size%2 == 0) {
        max_size /= 2;
        ++m_max_order;
        if (max_size <= max_block_size) {
            break;
        }
    }
    if (max_size > max_block_size) {
        amrex::Abort("DArena: incompatible max_size and max_block_size");
    }
    if (m_max_order > m_max_max_order) {
        amrex::Abort("DArena: too many orders");
    }
    //    const std::size_t a = alignof(std::max_align_t);
    const std::size_t a = 16;
    m_block_size = (max_size/a)*a;  // proper alignment
    m_max_size = m_block_size * (1u << m_max_order);

    if (amrex::Verbose()) {
        amrex::Print() << "DArena: Allocating " << m_max_size << " bytes\n";
    }

    m_baseptr = (char*)allocate_system(m_max_size);
    m_free[m_max_order].insert(0);
}

DArena::~DArena ()
{
    deallocate_system(m_baseptr);
}

void*
DArena::alloc (std::size_t nbytes)
{
    if (nbytes == 0) return nullptr; // behavior different from the standard

    // We need to allocate this many blocks
    unsigned int nblocks = (nbytes+m_block_size-1)/m_block_size;

    int order = 0;
    unsigned int ntmp = nblocks;
    while (ntmp >>= 1) {
        ++order;
    }
    if ((nblocks & (nblocks - 1)) != 0) { // nblocks is NOT power of 2
        ++order;
    }
    // We have found the miminal order that satisfies 2^order >= nblocks

    std::lock_guard<std::mutex> lock(m_mutex);

    std::ptrdiff_t offset = allocate_order(order);
    if (offset >= 0) {
        offset *= m_block_size; // # of order 0 blocks -> # of bytes
        m_used.insert({offset,order});
        return m_baseptr + offset;
    } else {
        if (amrex::Verbose()) {
            if (!warning_printed) {
                amrex::AllPrint() << "WARNING: DArena on proc. "
                                  << ParallelDescriptor::MyProc()
                                  << " calls system allocator to allocate "
                                  << nbytes << " bytes.\n"
                                  << "         DArena free memory size "
                                  << freeMem() << " bytes\n";
                warning_printed = true;
            }
        }
        return allocate_system(nbytes); // use the system malloc as backup.
    }
}

void
DArena::free (void* p)
{
    std::lock_guard<std::mutex> lock(m_mutex);

    std::ptrdiff_t offset = (char*)p - m_baseptr;
    auto r = m_used.find(offset);
    if (r != m_used.end()) {
        int order = r->second;
        m_used.erase(r);
        offset /= m_block_size;  // # of bytes -> # of order 0 blocks
        deallocate_order(order, offset);
    } else {
        deallocate_system(p);
    }
}

std::ptrdiff_t
DArena::allocate_order (int order)
{
    std::ptrdiff_t offset = -1;
    if (!m_free[order].empty()) {
        offset = *m_free[order].cbegin();
        m_free[order].erase(m_free[order].cbegin());
    } else if (order < m_max_order) {
        // We have to allocate from higher order that has larger blocks
        offset = allocate_order(order+1);
        if (offset >= 0) {
            // We have succeed. We need to split and put our buddy into free buckets
            m_free[order].insert(make_buddy(order,offset));
        }
    }

    return offset;
}

void
DArena::deallocate_order (int order, std::ptrdiff_t offset)
{
    std::ptrdiff_t buddy = make_buddy(order,offset);
    auto r = m_free[order].find(buddy);
    if (r != m_free[order].end()) {
        m_free[order].erase(r);
        deallocate_order(order+1,std::min(offset,buddy));
    } else {
        m_free[order].insert(offset);
    }
}

std::size_t
DArena::freeMem () const
{
    std::size_t r = 0;
    std::size_t bs = 1;
    for (int order = 0; order <= m_max_order; ++order) {
        r += bs*m_free[order].size();
        bs *= 2;
    }
    return r*m_block_size;
}

}
