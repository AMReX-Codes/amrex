#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <climits>

#include <AMReX_TagBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ccse-mpi.H>
#include <AMReX_iMultiFab.H>
#include <AMReX_MultiFabUtil.H>

namespace amrex {

TagBox::TagBox () noexcept {}

TagBox::TagBox (Arena* ar) noexcept
    : BaseFab<TagBox::TagType>(ar)
{}

TagBox::TagBox (const Box& bx, int n, Arena* ar)
    : BaseFab<TagBox::TagType>(bx,n,ar)
{}

TagBox::TagBox (const Box& bx, int n, bool alloc, bool shared, Arena* ar)
    : BaseFab<TagBox::TagType>(bx,n,alloc,shared,ar)
{}

TagBox::TagBox (const TagBox& rhs, MakeType make_type, int scomp, int ncomp)
    : BaseFab<TagBox::TagType>(rhs,make_type,scomp,ncomp)
{}

void
TagBox::coarsen (const IntVect& ratio, const Box& cbox) noexcept
{
    BL_ASSERT(nComp() == 1);
    Array4<char const> const& farr = this->const_array();

    TagBox cfab(cbox, 1, The_Arena());
    Elixir eli = cfab.elixir();
    Array4<char> const& carr = cfab.array();

    Box fdomain = domain;
    Dim3 r{1,1,1};
    AMREX_D_TERM(r.x = ratio[0];, r.y = ratio[1];, r.z = ratio[2]);

    AMREX_HOST_DEVICE_FOR_3D(cbox, i, j, k,
    {
        TagType t = TagBox::CLEAR;
        for (int koff = 0; koff < r.z; ++koff) {
            int kk = k*r.z + koff;
            for (int joff = 0; joff < r.y; ++joff) {
                int jj = j*r.y + joff;
                for (int ioff = 0; ioff < r.x; ++ioff) {
                    int ii = i*r.x + ioff;
                    if (fdomain.contains(IntVect(AMREX_D_DECL(ii,jj,kk)))) {
                        t = t || farr(ii,jj,kk);
                    }
                }
            }
        }
        carr(i,j,k) = t;
    });

#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        Gpu::dtod_memcpy_async(this->dataPtr(), cfab.dataPtr(), sizeof(TagType)*cbox.numPts());
    } else
#endif
    {
        std::memcpy(this->dataPtr(), cfab.dataPtr(), sizeof(TagType)*cbox.numPts());
    }
    this->domain = cbox;
}

void
TagBox::buffer (const IntVect& a_nbuff) noexcept
{
    Array4<char> const& a = this->array();
    Dim3 nbuf = a_nbuff.dim3();
    const auto lo = amrex::lbound(domain);
    const auto hi = amrex::ubound(domain);
    AMREX_HOST_DEVICE_FOR_3D(domain, i, j, k,
    {
        if (a(i,j,k) == TagBox::CLEAR) {
            bool to_buf = false;
            int imin = amrex::max(i-nbuf.x, lo.x);
            int jmin = amrex::max(j-nbuf.y, lo.y);
            int kmin = amrex::max(k-nbuf.z, lo.z);
            int imax = amrex::min(i+nbuf.x, hi.x);
            int jmax = amrex::min(j+nbuf.y, hi.y);
            int kmax = amrex::min(k+nbuf.z, hi.z);
            for (int kk = kmin; kk <= kmax && !to_buf; ++kk) {
            for (int jj = jmin; jj <= jmax && !to_buf; ++jj) {
            for (int ii = imin; ii <= imax && !to_buf; ++ii) {
                if (a(ii,jj,kk) == TagBox::SET) to_buf = true;
            }}}
            if (to_buf) a(i,j,k) = TagBox::BUF;
        }
    });
}

// DEPRECATED
Vector<int>
TagBox::tags () const noexcept
{
    Vector<int> ar(domain.numPts(), TagBox::CLEAR);

    const TagType* cptr = dataPtr();
    int*           iptr = ar.dataPtr();

    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        if (*cptr)
            *iptr = *cptr;
    }

    return ar;
}

// DEPRECATED
// Set values as specified by the array -- this only tags.
// It's an error if ar.length() != domain.numPts().
void
TagBox::tags (const Vector<int>& ar) noexcept
{
    BL_ASSERT(ar.size() == domain.numPts());

    TagType*   cptr = dataPtr();
    const int* iptr = ar.dataPtr();

    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        if (*iptr)
            *cptr = *iptr;
    }
}

// DEPRECATED
// Set values as specified by the array -- this tags and untags.
// It's an error if ar.length() != domain.numPts().
void
TagBox::tags_and_untags (const Vector<int>& ar) noexcept
{
    BL_ASSERT(ar.size() == domain.numPts());

    TagType*   cptr = dataPtr();
    const int* iptr = ar.dataPtr();

    // This clears as well as sets tags.
    for (int i = 0; i < ar.size(); i++, cptr++, iptr++)
    {
        *cptr = *iptr;
    }
}

// DEPRECATED
// Since a TagBox is a BaseFab<char>, we can use this utility
// function to allocate an integer array to have the same number 
// of elements as cells in tilebx
void 
TagBox::get_itags(Vector<int>& ar, const Box& tilebx) const noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Lbx[idim] = dlen[idim];
    }
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long Ntb = 1, stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Ltb[idim] = tilebx.length(idim);
        Ntb *= Ltb[idim];
        stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    if (ar.size() < Ntb) ar.resize(Ntb);
    
    const TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    int*                 iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
            const TagType* cptr = p0 + j*stride[1] + k*stride[2];
            for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
                if (*cptr) {
                    *iptr = *cptr;
                }
                else {
                    *iptr = TagBox::CLEAR;
                }
            }
        }
    }
}

// DEPRECATED
// Set values as specified by the array -- this only tags.
// only changes values in the tilebx region
void 
TagBox::tags (const Vector<int>& ar, const Box& tilebx) noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Lbx[idim] = dlen[idim];
    }
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Ltb[idim] = tilebx.length(idim);
        stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    const int*     iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
            TagType* cptr = p0 + j*stride[1] + k*stride[2];
            for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
                if (*iptr) *cptr = *iptr;
            }
        }
    }
}

// DEPRECATED
// Set values as specified by the array -- this tags and untags.
// only changes values in the tilebx region
void 
TagBox::tags_and_untags (const Vector<int>& ar, const Box& tilebx) noexcept
{
    auto dlen = length();
    int Lbx[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Lbx[idim] = dlen[idim];
    }
    
    Long stride[] = {1, Lbx[0], Long(Lbx[0])*Long(Lbx[1])};

    Long stb=0;
    int Ltb[] = {1,1,1};
    for (int idim=0; idim<AMREX_SPACEDIM; idim++) {
        Ltb[idim] = tilebx.length(idim);
        stb += stride[idim] * (tilebx.smallEnd(idim) - domain.smallEnd(idim));
    }
    
    TagType* const p0   = dataPtr() + stb;  // +stb to the lower corner of tilebox
    const int*     iptr = ar.dataPtr();

    for (int k=0; k<Ltb[2]; k++) {
        for (int j=0; j<Ltb[1]; j++) {
            TagType* cptr = p0 + j*stride[1] + k*stride[2];
            for (int i=0; i<Ltb[0]; i++, cptr++, iptr++) {
                *cptr = *iptr;
            }
        }
    }
}

TagBoxArray::TagBoxArray (const BoxArray& ba,
			  const DistributionMapping& dm,
                          int             _ngrow)
    :
    FabArray<TagBox>(ba,dm,1,_ngrow,MFInfo(),DefaultFabFactory<TagBox>())
{
    setVal(TagBox::CLEAR);
}

TagBoxArray::TagBoxArray (const BoxArray& ba,
			  const DistributionMapping& dm,
                          const IntVect&  _ngrow)
    :
    FabArray<TagBox>(ba,dm,1,_ngrow,MFInfo(),DefaultFabFactory<TagBox>())
{
    setVal(TagBox::CLEAR);
}

void 
TagBoxArray::buffer (const IntVect& nbuf)
{
    AMREX_ASSERT(nbuf.allLE(n_grow));

    if (nbuf.max() > 0)
    {
#ifdef _OPENMP
#pragma omp parallel if (Gpu::notInLaunchRegion())
#endif
       for (MFIter mfi(*this); mfi.isValid(); ++mfi) {
           get(mfi).buffer(nbuf);
       }
    }
}

void
TagBoxArray::mapPeriodicRemoveDuplicates (const Geometry& geom)
{
    BL_PROFILE("TagBoxArray::mapPRD");

    if (Gpu::inLaunchRegion())
    {
        // There is not atomicAdd for char.  So we have to use int.
        iMultiFab itag = amrex::cast<iMultiFab>(*this);
        iMultiFab tmp(boxArray(),DistributionMap(),1,nGrowVect());
        tmp.setVal(0);
        tmp.ParallelAdd(itag, 0, 0, 1, nGrowVect(), nGrowVect(), geom.periodicity());

        // We need to keep tags in periodic boundary
        const auto owner_mask = amrex::OwnerMask(tmp, Periodicity::NonPeriodic(), nGrowVect());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
            Box const& box = mfi.fabbox();
            Array4<TagType> const& tag =this->array(mfi);
            Array4<int const> const& tmptag = tmp.const_array(mfi);
            Array4<int const> const& msk = owner_mask->const_array(mfi);
            amrex::ParallelFor(box,
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
            {
                if (msk(i,j,k)) {
                    tag(i,j,k) = static_cast<char>(tmptag(i,j,k));
                } else {
                    tag(i,j,k) = TagBox::CLEAR;
                }
            });
        }
    }
    else
    {
        TagBoxArray tmp(boxArray(),DistributionMap(),nGrowVect()); // note that tmp is filled w/ CLEAR.
        tmp.ParallelAdd(*this, 0, 0, 1, nGrowVect(), nGrowVect(), geom.periodicity());

        // We need to keep tags in periodic boundary
        const auto owner_mask = amrex::OwnerMask(tmp, Periodicity::NonPeriodic(), nGrowVect());
#ifdef _OPENMP
#pragma omp parallel
#endif
        for (MFIter mfi(tmp); mfi.isValid(); ++mfi) {
            Box const& box = mfi.fabbox();
            Array4<TagType> const& tag = tmp.array(mfi);
            Array4<int const> const& msk = owner_mask->const_array(mfi);
            AMREX_LOOP_3D(box, i, j, k,
            {
                if (!msk(i,j,k)) tag(i,j,k) = TagBox::CLEAR;
            });
        }

        std::swap(*this, tmp);
    }
}

void
TagBoxArray::local_collate_cpu (Vector<IntVect>& v) const
{
    if (this->local_size() == 0) return;

    Vector<int> count(this->local_size());
#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Array4<char const> const& arr = this->const_array(fai);
        Box const& bx = fai.fabbox();
        int c = 0;
        AMREX_LOOP_3D(bx,i,j,k,
        {
            if (arr(i,j,k) != TagBox::CLEAR) ++c;
        });
        count[fai.LocalIndex()] = c;
    }

    Vector<int> offset(count.size()+1);
    offset[0] = 0;
    std::partial_sum(count.begin(), count.end(), offset.begin()+1);

    v.resize(offset.back());

    if (v.empty()) return;

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        int li = fai.LocalIndex();
        if (count[li] > 0) {
            IntVect* p = v.data() + offset[li];
            Array4<char const> const& arr = this->const_array(fai);
            Box const& bx = fai.fabbox();
            AMREX_LOOP_3D(bx,i,j,k,
            {
                if (arr(i,j,k) != TagBox::CLEAR) {
                    *p++ = IntVect(AMREX_D_DECL(i,j,k));
                }
            });
        }
    }
}

#ifdef AMREX_USE_GPU
void
TagBoxArray::local_collate_gpu (Vector<IntVect>& v) const
{
    const int nfabs = this->local_size();
    if (nfabs == 0) return;

    constexpr int block_size = 128;
    Vector<int> nblocks(nfabs);
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        Box const& bx = fai.fabbox();
        nblocks[fai.LocalIndex()] = (bx.numPts() + block_size-1) / block_size;
    }
    Vector<int> blockoffset(nblocks.size()+1);
    blockoffset[0] = 0;
    std::partial_sum(nblocks.begin(), nblocks.end(), blockoffset.begin()+1);
    int ntotblocks = blockoffset.back();

    PODVector<int,DeviceArenaAllocator<int> > dv_ntags(ntotblocks);

    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const int li = fai.LocalIndex();
        int* ntags = dv_ntags.data() + blockoffset[li];
        const int ncells = fai.fabbox().numPts();
        const char* tags = (*this)[fai].dataPtr();
#ifdef AMREX_USE_DPCPP
        amrex::launch(nblocks[li], block_size, sizeof(int)*Gpu::Device::warp_size,
                      Gpu::Device::gpuStream(),
        [=] AMREX_GPU_DEVICE (Gpu::Handler const& h) noexcept
        {
            int bid = h.item.get_group_linear_id();
            int tid = h.item.get_local_id(0);
            int icell = h.item.get_global_id(0);

            int t = 0;
            if (icell < ncells && tags[icell] != TagBox::CLEAR) {
                t = 1;
            }

            t = Gpu::blockReduce<Gpu::Device::warp_size>
                (t, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0, h);
            if (tid == 0) {
                ntags[bid] = t;
            }
        });
#else
        amrex::launch(nblocks[li], block_size, Gpu::Device::gpuStream(),
        [=] AMREX_GPU_DEVICE () noexcept
        {
            int bid = blockIdx.x;
            int tid = threadIdx.x;
            int icell = blockDim.x*blockIdx.x+threadIdx.x;

            int t = 0;
            if (icell < ncells && tags[icell] != TagBox::CLEAR) {
                t = 1;
            }

            t = Gpu::blockReduce<Gpu::Device::warp_size>
                (t, Gpu::warpReduce<Gpu::Device::warp_size,int,amrex::Plus<int> >(), 0);
            if (tid == 0) {
                ntags[bid] = t;
            }
        });
#endif
    }

    PODVector<int,PinnedArenaAllocator<int> > hv_ntags(ntotblocks);
    Gpu::dtoh_memcpy(hv_ntags.data(), dv_ntags.data(), ntotblocks*sizeof(int));

    PODVector<int,PinnedArenaAllocator<int> > hv_tags_offset(ntotblocks+1);
    hv_tags_offset[0] = 0;
    std::partial_sum(hv_ntags.begin(), hv_ntags.end(), hv_tags_offset.begin()+1);
    int ntotaltags = hv_tags_offset.back();

    if (ntotaltags == 0) return;

    PODVector<int,DeviceArenaAllocator<int> > dv_tags_offset(ntotblocks);
    int* dp_tags_offset = dv_tags_offset.data();
    Gpu::htod_memcpy(dp_tags_offset, hv_tags_offset.data(), ntotblocks*sizeof(int));
#ifdef AMREX_USE_DPCPP
    Gpu::synchronize();
#endif

    PODVector<IntVect,DeviceArenaAllocator<IntVect> > dv_tags(ntotaltags);
    IntVect* dp_tags = dv_tags.data();

    int iblock = 0;
    for (MFIter fai(*this); fai.isValid(); ++fai)
    {
        const int li = fai.LocalIndex();
        int iblock_begin = iblock;
        int iblock_end = iblock + nblocks[li];
        iblock = iblock_end;
        int count = 0;
        for (int ib = iblock_begin; ib < iblock_end; ++ib) {
            count += hv_ntags[ib];
        }
        if (count > 0) {
            Box const& bx = fai.fabbox();
            const auto lo  = amrex::lbound(bx);
            const auto len = amrex::length(bx);
            const int ncells = bx.numPts();
            const char* tags = (*this)[fai].dataPtr();
#ifdef AMREX_USE_DPCPP
            amrex::launch(nblocks[li], block_size, sizeof(unsigned int), Gpu::Device::gpuStream(),
            [=] AMREX_GPU_DEVICE (Gpu::Handler const& h) noexcept
            {
                int bid = h.item.get_group(0);
                int tid = h.item.get_local_id(0);
                int icell = h.item.get_global_id(0);

                unsigned int* shared_counter = (unsigned int*)h.local;
                if (tid == 0) {
                    *shared_counter = 0;
                }
                h.item.barrier(sycl::access::fence_space::local_space);

                if (icell < ncells && tags[icell] != TagBox::CLEAR) {
                    unsigned int itag = Gpu::Atomic::Inc<sycl::access::address_space::local_space>
                        (shared_counter, 20480u);
                    IntVect* p = dp_tags + dp_tags_offset[iblock_begin+bid];
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    p[itag] = IntVect(AMREX_D_DECL(i,j,k));
                }
            });
#else
            amrex::launch(nblocks[li], block_size, sizeof(unsigned int), Gpu::Device::gpuStream(),
            [=] AMREX_GPU_DEVICE () noexcept
            {
                int bid = blockIdx.x;
                int tid = threadIdx.x;
                int icell = blockDim.x*blockIdx.x+threadIdx.x;

                Gpu::SharedMemory<unsigned int> gsm;
                unsigned int * shared_counter = gsm.dataPtr();
                if (tid == 0) {
                    *shared_counter = 0;
                }
                __syncthreads();

                if (icell < ncells && tags[icell] != TagBox::CLEAR) {
                    unsigned int itag = Gpu::Atomic::Inc(shared_counter, blockDim.x);
                    IntVect* p = dp_tags + dp_tags_offset[iblock_begin+bid];
                    int k =  icell /   (len.x*len.y);
                    int j = (icell - k*(len.x*len.y)) /   len.x;
                    int i = (icell - k*(len.x*len.y)) - j*len.x;
                    i += lo.x;
                    j += lo.y;
                    k += lo.z;
                    p[itag] = IntVect(AMREX_D_DECL(i,j,k));
                }
            });
#endif
        }
    }

    v.resize(ntotaltags);
    Gpu::dtoh_memcpy(v.data(), dp_tags, ntotaltags*sizeof(IntVect));
}
#endif

void
TagBoxArray::collate (Vector<IntVect>& TheGlobalCollateSpace) const
{
    BL_PROFILE("TagBoxArray::collate()");

    Vector<IntVect> TheLocalCollateSpace;
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        local_collate_gpu(TheLocalCollateSpace);
    } else
#endif
    {
        local_collate_cpu(TheLocalCollateSpace);
    }

    Long count = TheLocalCollateSpace.size();

    //
    // The total number of tags system wide that must be collated.
    //
    Long numtags = count;
    ParallelDescriptor::ReduceLongSum(numtags);

    if (numtags == 0) {
        TheGlobalCollateSpace.clear();
        return;
    } else if (numtags > static_cast<Long>(std::numeric_limits<int>::max())) {
        // xxxxx todo
        amrex::Abort("TagBoxArray::collate: Too many tags. Using a larger blocking factor might help. Please file an issue on github");
    }

#ifdef BL_USE_MPI
    //
    // On I/O proc. this holds all tags after they've been gather'd.
    // On other procs. non-mempty signals size is not zero.
    //
    if (ParallelDescriptor::IOProcessor()) {
        TheGlobalCollateSpace.resize(numtags);
    } else {
        TheGlobalCollateSpace.resize(1);
    }

    //
    // Tell root CPU how many tags each CPU will be sending.
    //
    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    const std::vector<int>& countvec = ParallelDescriptor::Gather(static_cast<int>(count),
                                                                  IOProcNumber);
    std::vector<int> offset(countvec.size(),0);
    if (ParallelDescriptor::IOProcessor()) {
        for (int i = 1, N = offset.size(); i < N; i++) {
	    offset[i] = offset[i-1] + countvec[i-1];
	}
    }
    //
    // Gather all the tags to IOProcNumber into TheGlobalCollateSpace.
    //
    const IntVect* psend = (count > 0) ? TheLocalCollateSpace.data() : nullptr;
    IntVect* precv = TheGlobalCollateSpace.data();
    ParallelDescriptor::Gatherv(psend, count, precv, countvec, offset, IOProcNumber);

#else
    TheGlobalCollateSpace = std::move(TheLocalCollateSpace);
#endif
}

void
TagBoxArray::setVal (const BoxList& bl, TagBox::TagVal val)
{
    BoxArray ba(bl);
    setVal(ba,val);
}

void
TagBoxArray::setVal (const BoxDomain& bd, TagBox::TagVal val)
{
    setVal(bd.boxList(),val);
}

void
TagBoxArray::setVal (const BoxArray& ba, TagBox::TagVal val)
{
    Vector<Array4BoxTag<char> > tags;
    bool run_on_gpu = Gpu::inLaunchRegion();
#ifdef _OPENMP
#pragma omp parallel if (!run_on_gpu)
#endif
    {
        std::vector< std::pair<int,Box> > isects;
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            TagBox& fab = (*this)[mfi];
            Array4<char> const& arr = this->array(mfi);
            ba.intersections(mfi.fabbox(), isects);
            for (const auto& is : isects) {
                Box const& b = is.second;
                if (run_on_gpu) {
                    tags.push_back({arr,b});
                } else {
                   fab.setVal<RunOn::Host>(val,b);
                }
            }
        }
    }

#ifdef AMREX_USE_GPU
    amrex::ParallelFor(tags, 1,
    [=] AMREX_GPU_DEVICE (int i, int j, int k, int /*n*/, Array4<char> const& a) noexcept
    {
        a(i,j,k) = val;
    });
#endif
}

void
TagBoxArray::coarsen (const IntVect & ratio)
{
    // If team is used, all team workers need to go through all the fabs,
    // including ones they don't own.
    int teamsize = ParallelDescriptor::TeamSize();
    unsigned char flags = (teamsize == 1) ? 0 : MFIter::AllBoxes;

    IntVect new_n_grow;
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        new_n_grow[idim] = (n_grow[idim]+ratio[idim]-1)/ratio[idim];
    }

#if defined(_OPENMP)
#pragma omp parallel if (teamsize == 1 && Gpu::notInLaunchRegion())
#endif
    for (MFIter mfi(*this,flags); mfi.isValid(); ++mfi)
    {
        Box const& cbox = amrex::grow(amrex::coarsen(mfi.validbox(),ratio),new_n_grow);
        this->fabPtr(mfi)->coarsen(ratio,cbox);
    }

    boxarray.coarsen(ratio);
    n_grow = new_n_grow;
}

bool
TagBoxArray::hasTags (Box const& a_bx) const
{
    bool has_tags = false;
#ifdef AMREX_USE_GPU
    if (Gpu::inLaunchRegion()) {
        ReduceOps<ReduceOpLogicalOr> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;

        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            Box const& b = a_bx & mfi.fabbox();
            if (b.ok()) {
                const auto& arr = this->const_array(mfi);
                reduce_op.eval(b, reduce_data,
                [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
                {
                    int tr = arr(i,j,k) != TagBox::CLEAR;
                    return {tr};
                });
            }
        }

        ReduceTuple hv = reduce_data.value();
        has_tags = static_cast<bool>(amrex::get<0>(hv));
    } else
#endif
    {
#ifdef _OPENMP
#pragma omp parallel reduction(||:has_tags)
#endif
        for (MFIter mfi(*this); mfi.isValid(); ++mfi)
        {
            Box const& b = a_bx & mfi.fabbox();
            if (b.ok()) {
                Array4<char const> const& arr = this->const_array(mfi);
                AMREX_LOOP_3D(b, i, j, k,
                {
                    has_tags = has_tags || (arr(i,j,k) != TagBox::CLEAR);
                });
            }
        }
    }

    ParallelAllReduce::Or(has_tags, ParallelContext::CommunicatorSub());
    return has_tags;
}

}

