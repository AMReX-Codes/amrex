//BL_COPYRIGHT_NOTICE

//
// $Id: CArena.cpp,v 1.4 1998-02-07 23:41:57 lijewski Exp $
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


#include <ParallelDescriptor.H>

static CArena The_Static_FAB_CArena;

extern CArena* The_FAB_CArena = &The_Static_FAB_CArena;

//
// These routines implement a COMPACTING version of Kernighan
// and Ritchie's memory manager.  Chunks of memory are allocated
// by the system memory manager.  The union of all these chunks
// are broken up into a list of blocks, each of which is either
// in use by the calling program or free to be used.
// 
// If the block is in use it contains an extra header at the
// front of the block that records the size of the used block and
// pointers to the next and previous blocks.  It also contains
// the address of the pointer in the user program that points
// to this used space (for use in compacting).
//
// If the block is free, the first header in the block
// contains the size of the free block, pointers to the
// next and previous free blocks (the free list) and
// pointers to the next and previous blocks (used or free).
//
// CKR_INIT:  sets the size of the default memory chunks.
//
// CKR_MALLOC:  searches the free list for a block of memory
// large enough to satisfy the request.  It modifies the block
// list and free list to indicate that this space is taken.
// If there is no free block large enough to satisfy the
// request then the memory is compacted (creating larger free
// blocks).  If this still fails to satisfy the request then
// a new chunk is allocated from the system.  The size of this
// new chunk is the maximum of the size of the requested
// memory and the default size of the chunk.  If the system
// cannot provide a chunk that size, the memory manager fails
// and execution is halted.
//
// CKR_FREE: returns space to the free list.  The header
// in front of the returned memory address contains the
// information about that block.  If the previous or next
// blocks are adjacent to this block and are free, it will
// compact the blocks into one large free block.
//
// variable which affects the compaction policy:
//
// c_hunk_size sets the size of memory hunks which morecore() will try to
// allocate.  It is measured in units of sizeof(Header)
//

bool CArena::m_verbose = false;

CArena::CArena (size_t size)
{
    base.s.size = 0;
    base.s.frprev = base.s.frnext = freep = &base;
    base.s.blnext = base.s.blprev = &base;
    highest_mem = lowest_mem = &base;

    if (size == 0)
    {
#if   defined(BOXLIB_ARCH_IEEE)
        c_hunk_size = (8000000/sizeof(Header));
#elif defined(BOXLIB_ARCH_CRAY)
        c_hunk_size = (40000000/sizeof(Header));
#endif
    }
    else
    {
        c_hunk_size = size/sizeof(Header);
    }
}

CArena::~CArena ()
{
    for (int i = 0; i < psaved.length(); ++i)
    {
        ::operator delete(psaved[i]);
    }
}

//
// Static statistic recording variables.
//
static size_t bytes_alloc;
static size_t bytes_free;
static int    calls_malloc;   // Actually to `::operator new()'
static int    calls_free;     // Actually to `::operator delete()'
static int    calls_compact;
static int    calls_more;

void*
CArena::alloc (size_t nbytes,
               void** ind)
{
    assert(!(ind == 0));

    calls_malloc++;
    unsigned nunits = (nbytes+sizeof(Header)-1)/ sizeof(Header)+1;
    Header* prevp = freep;
    //
    // Search the free list.
    //
    for (Header* p = prevp->s.frnext; ; prevp = p, p = p->s.frnext)
    {
        if (p > highest_mem || p < lowest_mem)
            BoxLib::Error("memory corrupted");

        if (p->s.size >= nunits)
        {
            if (p->s.size == nunits)
            {
                //
                // Exactly.
                //
                prevp->s.frnext = p->s.frnext;
                p->s.frnext->s.frprev = prevp;
            }
            else
            {
                //
                // Not exact fit, must make new block.
                // freepart is the hi part, busy part is lo.
                //
                Header* freepart = p + nunits;
                freepart->s.frnext = p->s.frnext;
                p->s.frnext->s.frprev = freepart;
                freepart->s.frprev = p->s.frprev;
                p->s.frprev->s.frnext = freepart;
                freepart->s.size = p->s.size - nunits;
                freepart->s.blnext = p->s.blnext;
                freepart->s.blprev = p;
                //
                // Also correct pointer of next block.
                //
                freepart->s.blnext->s.blprev = freepart;
                p->s.size = nunits;
                p->s.blnext = freepart;
            }
            p->s.frnext = 0;    // Indicates is not free.
            p->s.frprev = 0;
            freep = prevp;
            *ind = p+1;
            p->s.callback = ind;
            bytes_free -= p->s.size*sizeof(Header);
            return *ind;
        }

        if (p == freep && (p = morecore(nunits)) == 0)
        {
            BoxLib::Error("Out Of Memory");
            *ind = 0;
            return *ind;
        }
    }
}

CArena::Header*
CArena::morecore (size_t nu)
{
    calls_more++;

    if (nu < c_hunk_size)
        nu = c_hunk_size;

    unsigned int kk = nu*sizeof(Header);
    void* vp = ::operator new(kk);
    //
    // Now save my pointer for later deletion.
    //
    int olen = psaved.length();
    psaved.resize(olen+1);
    psaved[olen] = vp;
    bytes_alloc += kk;
    Header* up = (Header*) vp;
    if (up < lowest_mem)
        lowest_mem = up;
    if (&up[nu-1] > highest_mem)
        highest_mem = &up[nu-1];
    up->s.size = kk/sizeof(Header);
    //
    // Install on end of block list.
    //
    up->s.blnext = &base;
    up->s.blprev = base.s.blprev;
    base.s.blprev->s.blnext = up;
    base.s.blprev = up;
    CArena::free(up+1);
    if (m_verbose)
    {
        clog << "CArena::morecore(): new = "
             << bytes_alloc
             << " free = "
             << bytes_free
             << " ("
             << 100.0*double(bytes_free)/double(bytes_alloc)
             << " % of total)"
             << '\n';
    }
    return freep;
}

void
CArena::free (void* ap)
{
    if (ap == 0)
        return;

    calls_free++;
    Header *bp = (Header*) ap, *p = freep;
    bp -= 1; // Point to block header.
    bytes_free += bp->s.size*sizeof(Header);
    //
    // Search the free list.
    //
    for ( ; !(bp > p && bp < p->s.frnext); p = p->s.frnext)
    {
        if (p > highest_mem || p < lowest_mem)
            BoxLib::Error("memory corrupted");

        if (p >= p->s.frnext && (bp > p || bp < p->s.frnext))
            //
            // Freed block at start or end.
            //
            break;
    }
    if (bp + bp->s.size == p->s.frnext)
    {
        //
        // Join to upper nbr.
        //
        Header* above = p->s.frnext;
        bp->s.size += above->s.size;
        bp->s.frnext = above->s.frnext;
        bp->s.frnext->s.frprev = bp;
        bp->s.blnext = above->s.blnext;
        (bp->s.blnext)->s.blprev = bp;
    }
    else
    {
        bp->s.frnext = p->s.frnext;
        bp->s.frnext->s.frprev = bp;
    }
    if (p + p->s.size == bp)
    {
        //
        // Join to lower nbr.
        //
        p->s.size += bp->s.size;
        p->s.frnext = bp->s.frnext;
        p->s.frnext->s.frprev = p;
        p->s.blnext = bp->s.blnext;
        (bp->s.blnext)->s.blprev = p;
    }
    else
    {
        p->s.frnext = bp;
        bp->s.frprev = p;
    }
    freep = p;
}


//
// Returns pointer to next block to work on.  If blocks coalesce, then
// points at lower free block.  If blocks do not coalesce because they
// are not contiguous, then points at higher block.
//

CArena::Header*
CArena::freeFree (Header* lower,
                  Header* upper)
{
    //
    // Check that are really free.
    //
    if (upper->s.frnext == 0 || upper->s.frprev == 0 ||
        lower->s.frnext == 0 || lower->s.frprev == 0)
    {
        BoxLib::Error("free block is actually used");
    }
    //
    // Check that are adjacent.
    //
    if (upper != (lower + lower->s.size))
        return upper;
    //
    // Reset values in lower.
    //
    lower->s.size  += upper->s.size;
    lower->s.frnext = upper->s.frnext;
    lower->s.blnext = upper->s.blnext;
    //
    // Correct pointers that pointed to upper.
    //
    assert(upper->s.frnext->s.frprev == upper);
    assert(upper->s.blnext->s.blprev == upper);
    upper->s.frnext->s.frprev = lower;
    upper->s.blnext->s.blprev = lower;

    return lower;
}

//
// Returns pointer to next block to work on.  If blocks are contiguous, so
// can do an exchange, will point at (new location of) free block.  If
// blocks are not contiguous, will point at next free block.
//

CArena::Header*
CArena::freeUsed (Header* freeblock,
                  Header* busyblock)
{
    if (busyblock->s.frnext || busyblock->s.frprev)
        BoxLib::Error("busy block is actually free");

    if (busyblock != (freeblock + freeblock->s.size))
        //
        // Are not contiguous.
        //
        return freeblock->s.frnext;

    Header oldfree = *freeblock;    // copy stuff
    Header oldbusy = *busyblock;
    //
    // Save pointer to busy memory.
    //
    void* busymem = busyblock+1;
    //
    // Reset block pointers.
    //
    busyblock = freeblock;
    freeblock = busyblock + oldbusy.s.size;
    //
    // New values for the lower block, which will be busy.
    //
    busyblock->s.frnext = busyblock->s.frprev = 0;
    busyblock->s.size   = oldbusy.s.size;
    busyblock->s.blnext = freeblock;
    busyblock->s.blprev = oldfree.s.blprev;
    //
    // Change the callback pointer.
    //
    busyblock->s.callback    = oldbusy.s.callback;
    *(busyblock->s.callback) = busyblock+1;
    //
    // Move the memory.
    //
    memmove(busyblock+1, busymem, (oldbusy.s.size-1)*sizeof(Header));
    //
    // New values for the upper block, which will be free.
    //
    freeblock->s.frnext = oldfree.s.frnext;
    freeblock->s.frprev = oldfree.s.frprev;
    freeblock->s.size   = oldfree.s.size;
    freeblock->s.blnext = oldbusy.s.blnext;
    freeblock->s.blprev = busyblock;
    //
    // Reset free neighbors.
    //
    freeblock->s.frnext->s.frprev = freeblock;
    freeblock->s.frprev->s.frnext = freeblock;
    //
    // Reset block neighbors.
    //
    freeblock->s.blnext->s.blprev = freeblock;

    return freeblock;
}

void
CArena::compact ()
{
    calls_compact++;

    assert(CArena::ok(0));

    for (Header* p = base.s.frnext; p != &base ; )
    {
        assert(CArena::ok(0));
        Header* nxtblk = p->s.blnext;
        if (nxtblk == &base)
            break;
        if (nxtblk->s.frnext)
        {
            p = freeFree(p,nxtblk);
        }
        else
        {
            p = freeUsed(p,nxtblk);
        }
        assert(CArena::ok(0));
    }
    //
    // freep must be reset so it points to an element of free list.
    //
    freep = &base;
}

void
CArena::dump (const char* what,
              void*       p)
{
    ofstream dump("mem_map", ios::out);

    dump << "CArena::ok(): "
         << what
         << " chain corrupted near "
         << p
         << '\n';

    stat(dump,3);
    dump.close();
    BoxLib::Abort();
}

//
// level = 0 means do not write out memory map.
//

bool
CArena::ok (int level)
{
    for (Header* p = base.s.blprev; p != &base; p = p->s.blprev)
    {
        Header* hnext = p->s.blnext;
        Header* hprev = p->s.blprev;

        if (hnext > highest_mem  ||
            hnext < lowest_mem   ||
            hprev > highest_mem  ||
            hprev < lowest_mem   ||
            hnext->s.blprev != p || 
            hprev->s.blnext != p)
        {
            cerr << "Block chain corrupted near " << p << '\n';
            if (level)
                dump("block", p);
            return false;
        }
    }
    for (Header* p = base.s.frprev; p != &base; p = p->s.frprev)
    {
        Header* hnext = p->s.frnext;
        Header* hprev = p->s.frprev;
        
        if (hnext > highest_mem  ||
            hnext < lowest_mem   ||
            hprev > highest_mem  ||
            hprev < lowest_mem   ||
            hnext->s.frprev != p || 
            hprev->s.frnext != p)
        {
            cerr << "Free chain corrupted near " << p << '\n';
            if (level)
                dump("free", p);
            return false;
        }
    }
    return true;
}

//
// level is bit encoded
// level&1 means turn on statistics
// level&2 means dump map in ascii
//

void
CArena::stat (ostream& os,
              int      level)
{
    Header* p;

    if (level&1)
    {
        //
        // Traverse block list to get stats.
        //
        os << calls_malloc  << " calls to malloc\n";
        os << calls_free    << " calls to free\n";
        os << calls_more    << " calls to morecore\n";
        os << calls_compact << " calls to compact\n";
        os << bytes_alloc   << " bytes allocated\n";
        os << bytes_free << " free bytes = "
           << (bytes_alloc==0 ? 0 : 100*double(bytes_free)/double(bytes_alloc))
           << " % of total" << '\n';
    }
    if (level&2)
    {
        //
        // Dump out ascii map.
        //
        os <<
            "    address"
            "       size"
            "  free-next"
            "  free-prev"
            " block-next"
            "   blk-prev"
            "   callback"
           << '\n';
    
        for (p = base.s.blnext; p != &base; p = p->s.blnext)
        {
            if (p > highest_mem || p < lowest_mem)
            {
                os << "pointer "
                   << p
                   << "is illegal, stopping dump"
                   << '\n';
                BoxLib::Abort();
            }
            os << setw(11) << p;
            os << setw(11) << p->s.size;
            os << setw(11) << p->s.frnext;
            os << setw(11) << p->s.frprev;
            os << setw(11) << p->s.blnext;
            os << setw(11) << p->s.blprev;
            os << setw(11) << p->s.callback;
            os << '\n';
        }
    }
}
