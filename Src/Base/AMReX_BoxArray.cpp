
#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>
#include <AMReX_MFIter.H>
#include <AMReX_BaseFab.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace amrex {

#ifdef BL_MEM_PROFILING
int  BARef::numboxarrays         = 0;
int  BARef::numboxarrays_hwm     = 0;
long BARef::total_box_bytes      = 0L;
long BARef::total_box_bytes_hwm  = 0L;
long BARef::total_hash_bytes     = 0L;
long BARef::total_hash_bytes_hwm = 0L;
#endif

bool    BARef::initialized = false;
bool BoxArray::initialized = false;

namespace {
    const int bl_ignore_max = 100000;
}

BARef::BARef () 
{ 
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif	    
}

BARef::BARef (size_t size) 
    : m_abox(size) 
{ 
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif	    
}
 
BARef::BARef (const Box& b)
{ 
    define(b); 
}

BARef::BARef (const BoxList& bl)
    : m_abox(bl.data())
{ 
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif	    
}

BARef::BARef (BoxList&& bl) noexcept
    : m_abox(std::move(bl.data()))
{ 
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif	    
}

BARef::BARef (std::istream& is)
{ 
    int ndims;
    define(is, ndims); 
}

BARef::BARef (const BARef& rhs) 
    : m_abox(rhs.m_abox) // don't copy hash
{
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif	    
}

BARef::~BARef ()
{
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
    updateMemoryUsage_hash(-1);
#endif	    
}

void
BARef::define (std::istream& is, int& ndims)
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    BL_ASSERT(m_abox.size() == 0);
    int           maxbox;
    unsigned long tmphash;
    is.ignore(bl_ignore_max, '(') >> maxbox >> tmphash;
    resize(maxbox);
    auto pos = is.tellg();
    {
        ndims = AMREX_SPACEDIM;
        char c1, c2;
        int itmp;
        is >> std::ws >> c1 >> std::ws >> c2;
        if (c1 == '(' && c2 == '(') {
            is >> itmp;
            ndims = 1;
#if (AMREX_SPACEDIM >= 2)
            is >> std::ws;
            int ic = is.peek();
            if (ic == static_cast<int>(',')) {
                is.ignore(BL_IGNORE_MAX, ',');
                is >> itmp;
                ++ndims;
#if (AMREX_SPACEDIM == 3)
                is >> std::ws;
                ic = is.peek();
                if (ic == static_cast<int>(',')) {
                    ++ndims;
                }
#endif
            }
#endif
        }
    }
    is.seekg(pos, std::ios_base::beg);
    for (Vector<Box>::iterator it = m_abox.begin(), End = m_abox.end(); it != End; ++it)
        is >> *it;
    is.ignore(bl_ignore_max, ')');
    if (is.fail())
        amrex::Error("BoxArray::define(istream&) failed");
}

void
BARef::define (const Box& bx)
{
    BL_ASSERT(m_abox.size() == 0);
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
#endif
    m_abox.push_back(bx);
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif
}

void
BARef::define (const BoxList& bl)
{
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
#endif
    m_abox = bl.data();
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif
}

void
BARef::define (BoxList&& bl) noexcept
{
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
#endif
    m_abox = std::move(bl.data());
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif
}

void 
BARef::resize (long n) {
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
    updateMemoryUsage_hash(-1);
#endif
    m_abox.resize(n);
    hash.clear();
    has_hashmap = false;
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(1);
#endif
}

#ifdef BL_MEM_PROFILING
void
BARef::updateMemoryUsage_box (int s)
{
    if (m_abox.size() > 1) {
	long b = amrex::bytesOf(m_abox);
	if (s > 0) {
	    total_box_bytes += b;
	    total_box_bytes_hwm = std::max(total_box_bytes_hwm, total_box_bytes);
	    ++numboxarrays;
	    numboxarrays_hwm = std::max(numboxarrays_hwm, numboxarrays);
	} else {
	    total_box_bytes -= b;
	    --numboxarrays;
	}
    }
}

void
BARef::updateMemoryUsage_hash (int s)
{
    if (hash.size() > 0) {
	long b = sizeof(hash);
	for (const auto& x: hash) {
	    b += amrex::gcc_map_node_extra_bytes
		+ sizeof(IntVect) + amrex::bytesOf(x.second);
	}
	if (s > 0) {
	    total_hash_bytes += b;
	    total_hash_bytes_hwm = std::max(total_hash_bytes_hwm, total_hash_bytes);
	} else {
	    total_hash_bytes -= b;
	}
    }
}
#endif

void
BARef::Initialize ()
{
    if (!initialized) {
	initialized = true;
#ifdef BL_MEM_PROFILING
	MemProfiler::add("BoxArray", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {total_box_bytes, total_box_bytes_hwm};
			 }));
	MemProfiler::add("BoxArrayHash", std::function<MemProfiler::MemInfo()>
			 ([] () -> MemProfiler::MemInfo {
			     return {total_hash_bytes, total_hash_bytes_hwm};
			 }));
	MemProfiler::add("BoxArray Innard", std::function<MemProfiler::NBuildsInfo()>
			 ([] () -> MemProfiler::NBuildsInfo {
			     return {numboxarrays, numboxarrays_hwm};
			 }));
#endif
    }

    amrex::ExecOnFinalize(BARef::Finalize);
}

void
BARef::Finalize ()
{
    initialized = false;
}

void
BoxArray::Initialize ()
{
    if (!initialized) {
	initialized = true;
	BARef::Initialize();
    }

    amrex::ExecOnFinalize(BoxArray::Finalize);
}

void
BoxArray::Finalize ()
{
    initialized = false;
}

BoxArray::BoxArray ()
    :
    m_transformer(new DefaultBATransformer()),
    m_typ(),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>())
{}

BoxArray::BoxArray (const Box& bx)
    :
    m_transformer(new DefaultBATransformer(bx.ixType())),
    m_typ(bx.ixType()),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>(amrex::enclosedCells(bx)))
{}

BoxArray::BoxArray (const BoxList& bl)
    :
    m_transformer(new DefaultBATransformer(bl.ixType())),
    m_typ(bl.ixType()),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>(bl))
{
    type_update();
}

BoxArray::BoxArray (BoxList&& bl) noexcept
    :
    m_transformer(new DefaultBATransformer(bl.ixType())),
    m_typ(bl.ixType()),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>(std::move(bl)))
{
    type_update();
}

BoxArray::BoxArray (size_t n)
    :
    m_transformer(new DefaultBATransformer()),
    m_typ(),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>(n))
{}

BoxArray::BoxArray (const Box* bxvec,
                    int        nbox)
    :
    m_transformer(new DefaultBATransformer(bxvec->ixType())),
    m_typ(bxvec->ixType()),
    m_crse_ratio(IntVect::TheUnitVector()),
    m_simple(true),
    m_ref(std::make_shared<BARef>(nbox))
{
    for (int i = 0; i < nbox; i++) {
        m_ref->m_abox[i] = amrex::enclosedCells(*bxvec++);
    }
}

BoxArray::BoxArray (const BoxArray& rhs, const BATransformer& trans)
    :
    m_transformer(trans.clone()),
    m_typ(trans.ixType()),
    m_crse_ratio(rhs.m_crse_ratio * trans.crseRatio()),
    m_simple(trans.simple()),
    m_ref(rhs.m_ref)
{
    BL_ASSERT(rhs.m_typ.cellCentered());  // rhs must be cell-centered.
    BL_ASSERT(rhs.m_simple);              // rhs must be simple!
    // But it's OK to have non-unit crse ratio. In that case we need to update the transformer
    m_transformer->setCrseRatio(m_crse_ratio);
}

BoxArray::BoxArray (const BoxArray& rhs)
    :
    m_transformer(rhs.m_transformer->clone()),
    m_typ(rhs.m_typ),
    m_crse_ratio(rhs.m_crse_ratio),
    m_simple(rhs.m_simple),
    m_ref(rhs.m_ref)
{}

BoxArray&
BoxArray::operator= (const BoxArray& rhs)
{
    m_transformer.reset(rhs.m_transformer->clone());
    m_typ = rhs.m_typ;
    m_crse_ratio = rhs.m_crse_ratio;
    m_simple = rhs.m_simple;
    m_ref = rhs.m_ref;
    return *this;
}

void
BoxArray::define (const Box& bx)
{
    clear();
    m_transformer->setIxType(bx.ixType());
    m_typ = bx.ixType();
    m_crse_ratio = IntVect::TheUnitVector();
    m_simple = true;
    m_ref->define(amrex::enclosedCells(bx));
}

void
BoxArray::define (const BoxList& bl)
{
    clear();
    m_ref->define(bl);
    m_crse_ratio = IntVect::TheUnitVector();
    m_simple = true;
    type_update();
}

void
BoxArray::define (BoxList&& bl) noexcept
{
    clear();
    m_ref->define(std::move(bl));
    m_crse_ratio = IntVect::TheUnitVector();
    m_simple = true;
    type_update();
}

void
BoxArray::clear ()
{
    m_transformer.reset(new DefaultBATransformer());
    m_ref.reset(new BARef());
}

void
BoxArray::resize (long len)
{
    uniqify();
    m_ref->resize(len);
}

long
BoxArray::numPts () const noexcept
{
    long result = 0;
    const int N = size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
    for (int i = 0; i < N; ++i)
    {
        result += (*this)[i].numPts();
    }
    return result;
}

double
BoxArray::d_numPts () const noexcept
{
    double result = 0;
    const int N = size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
    for (int i = 0; i < N; ++i)
    {
        result += (*this)[i].d_numPts();
    }
    return result;
}

int
BoxArray::readFrom (std::istream& is)
{
    BL_ASSERT(size() == 0);
    clear();
    m_simple = true;
    m_crse_ratio = IntVect::TheUnitVector();
    int ndims;
    m_ref->define(is, ndims);
    type_update();
    return ndims;
}

std::ostream&
BoxArray::writeOn (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << '(' << size() << ' ' << 0 << '\n';

    const int N = size();
    for (int i = 0; i < N; ++i) {
        os << (*this)[i] << '\n';
    }

    os << ')';

    if (os.fail())
        amrex::Error("BoxArray::writeOn(ostream&) failed");

    return os;
}

bool
BoxArray::operator== (const BoxArray& rhs) const noexcept
{
    if (m_simple && rhs.m_simple) {
        return m_typ == rhs.m_typ && m_crse_ratio == rhs.m_crse_ratio &&
            (m_ref == rhs.m_ref || m_ref->m_abox == rhs.m_ref->m_abox);
    } else {
        return m_simple == rhs.m_simple
            && m_typ == rhs.m_typ
            && m_crse_ratio == rhs.m_crse_ratio
            && m_transformer->equal(*rhs.m_transformer)
            && (m_ref == rhs.m_ref || m_ref->m_abox == rhs.m_ref->m_abox);
    }
}

bool
BoxArray::operator!= (const BoxArray& rhs) const noexcept
{
    return !operator==(rhs);
}

bool
BoxArray::operator== (const Vector<Box>& bv) const noexcept
{
    if (size() != bv.size()) return false;
    for (long i = 0; i < size(); ++i) {
        if (this->operator[](i) != bv[i]) return false;
    }
    return true;
}

bool
BoxArray::operator!= (const Vector<Box>& bv) const noexcept
{
    return !operator==(bv);
}

bool
BoxArray::CellEqual (const BoxArray& rhs) const noexcept
{
    return m_crse_ratio == rhs.m_crse_ratio
        && (m_ref == rhs.m_ref || m_ref->m_abox == rhs.m_ref->m_abox);
}

BoxArray&
BoxArray::maxSize (int block_size)
{
    return maxSize(IntVect(AMREX_D_DECL(block_size,block_size,block_size)));
}

BoxArray&
BoxArray::maxSize (const IntVect& block_size)
{
    if (!m_simple || m_crse_ratio != IntVect::TheUnitVector()) {
        uniqify();
    }
    BoxList blst(*this);
    blst.maxSize(block_size);
    const int N = blst.size();
    if (size() != N) { // If size doesn't change, do nothing.
        define(std::move(blst));
    }
    return *this;
}

BoxArray&
BoxArray::refine (int refinement_ratio)
{
    return refine(IntVect(AMREX_D_DECL(refinement_ratio,refinement_ratio,refinement_ratio)));
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
	BL_ASSERT(m_ref->m_abox[i].ok());
        m_ref->m_abox[i].refine(iv);
    }
    return *this;
}

bool
BoxArray::coarsenable(int refinement_ratio, int min_width) const
{
    return coarsenable(IntVect{refinement_ratio}, min_width);
}

bool
BoxArray::coarsenable(const IntVect& refinement_ratio, int min_width) const
{
    const long sz = size();
    if(size() == 0) return false;

    const Box& first = (*this)[0];
    bool res = first.coarsenable(refinement_ratio,min_width);
    if (res == false) return false;

#ifdef _OPENMP
#pragma omp parallel for reduction(&&:res)
#endif
    for (long ibox = 0; ibox < sz; ++ibox)
    {
        const Box& thisbox = (*this)[ibox];        
        res = res && thisbox.coarsenable(refinement_ratio,min_width);
    }        

    return res;
}

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    return coarsen(IntVect(AMREX_D_DECL(refinement_ratio,refinement_ratio,refinement_ratio)));
}

BoxArray&
BoxArray::coarsen (const IntVect& iv)
{
    m_crse_ratio *= iv;
    m_transformer->setCrseRatio(m_crse_ratio);
    return *this;
}

BoxArray&
BoxArray::growcoarsen (int n, const IntVect& iv)
{
    return growcoarsen(IntVect(n), iv);
}

BoxArray&
BoxArray::growcoarsen (IntVect const& ngrow, const IntVect& iv)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].grow(ngrow).coarsen(iv);
    }
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].grow(n);
    }
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].grow(iv);
    }
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].grow(dir, n_cell);
    }
    return *this;
}

BoxArray&
BoxArray::growLo (int dir,
                  int n_cell)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].growLo(dir, n_cell);
    }
    return *this;
}

BoxArray&
BoxArray::growHi (int dir,
                  int n_cell)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].growHi(dir, n_cell);
    }
    return *this;
}

BoxArray&
BoxArray::surroundingNodes ()
{
    
    return this->convert(IndexType::TheNodeType());
}

BoxArray&
BoxArray::surroundingNodes (int dir)
{
    IndexType typ = ixType();
    typ.setType(dir, IndexType::NODE);
    return this->convert(typ);
}

BoxArray&
BoxArray::enclosedCells ()
{
    return this->convert(IndexType::TheCellType());
}

BoxArray&
BoxArray::enclosedCells (int dir)
{
    IndexType typ = ixType();
    typ.setType(dir, IndexType::CELL);
    return this->convert(typ);
}

BoxArray&
BoxArray::convert (IndexType typ)
{
    m_transformer->setIxType(typ);
    m_typ = typ;
    return *this;
}

BoxArray&
BoxArray::convert (const IntVect& iv)
{
    IndexType typ(iv);
    m_transformer->setIxType(typ);
    m_typ = typ;
    return *this;
}

BoxArray&
BoxArray::convert (Box (*fp)(const Box&))
{
    BL_ASSERT(!(fp == 0));

    const int N = size();
    if (N > 0) {
        uniqify();

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; ++i) {
            set(i,fp((*this)[i]));
        }
    }
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].shift(dir, nzones);
    }
    return *this;
}

BoxArray&
BoxArray::shift (const IntVect& iv)
{
    uniqify();

    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++) {
        m_ref->m_abox[i].shift(iv);
    }
    return *this;
}

void
BoxArray::set (int        i,
               const Box& ibox)
{
    BL_ASSERT(m_simple && m_crse_ratio == IntVect::TheUnitVector());
    if (i == 0) {
        m_typ = ibox.ixType();
        m_transformer->setIxType(m_typ);
    }
    m_ref->m_abox[i] = amrex::enclosedCells(ibox);
}

Box
BoxArray::operator[] (const MFIter& mfi) const noexcept
{
    return (*this)[mfi.index()];
}

bool
BoxArray::ok () const
{
    const int N = size();
    if (N > 0)
    {
        for (int i = 0; i < N; ++i)
        {
            if(!(*this)[i].ok()) {
                return false;
            }
        }
    }
    return true;
}

bool
BoxArray::isDisjoint () const
{
    BL_ASSERT(ixType().cellCentered());

    std::vector< std::pair<int,Box> > isects;

    const int N = size();
    for (int i = 0; i < N; ++i)
    {
	intersections((*this)[i],isects);
        if ( isects.size() > 1 ) {
            return false;
        }
    }

    return true;
}

BoxList
BoxArray::boxList () const
{
    const int N = size();
    BoxList newb;
    newb.data().reserve(N);
    if (N > 0) {
	newb.set(ixType());
	for (int i = 0; i < N; ++i) {
	    newb.push_back((*this)[i]);
        }
    }
    return newb;
}

bool
BoxArray::contains (const IntVect& iv) const
{
    if (size() > 0) {
	return intersects(Box(iv,iv,ixType()));
    } else {
	return false;
    }
}

bool
BoxArray::contains (const Box& b, bool assume_disjoint_ba) const
{
    bool result = false;

    if (size() > 0)
    {
        BL_ASSERT(ixType() == b.ixType());

        std::vector< std::pair<int,Box> > isects;

        intersections(b,isects);

        if (isects.size() > 0)
        {
	    if (assume_disjoint_ba) {
                BL_ASSERT(isDisjoint());
		long nbx = b.numPts(), nisects = 0L;
		for (int i = 0, N = isects.size(); i < N; i++) {
		    nisects += isects[i].second.numPts();
		}
		result = nbx == nisects;
	    } else {
                Vector<char> vflag(b.numPts(), 1);
                BaseFab<char> fabflag(b, 1, vflag.data());
                for (int i = 0, N = isects.size(); i < N; i++) {
                    fabflag.setVal(0, isects[i].second, 0, 1);
                }
                for (const auto& x : vflag) {
                    if (x == 1) return false;
                }
                result = true;
	    }
        }
    }

    return result;
}

bool
BoxArray::contains (const BoxArray& ba, bool assume_disjoint_ba) const
{
    if (size() == 0) return false;

    if (!minimalBox().contains(ba.minimalBox())) return false;

    for (int i = 0, N = ba.size(); i < N; ++i) {
        if (!contains(ba[i],assume_disjoint_ba)) {
            return false;
	}
    }

    return true;
}

Box
BoxArray::minimalBox () const
{
    BL_ASSERT(m_simple);
    Box minbox;
    const int N = size();
    if (N > 0)
    {
#ifdef _OPENMP
	bool use_single_thread = omp_in_parallel();
	const int nthreads = use_single_thread ? 1 : omp_get_max_threads();
#else
	bool use_single_thread = true;
	const int nthreads = 1;
#endif
	if (use_single_thread)
	{
	    minbox = m_ref->m_abox[0];
	    for (int i = 1; i < N; ++i) {
		minbox.minBox(m_ref->m_abox[i]);
	    }
	}
	else
	{
	    Vector<Box> bxs(nthreads, m_ref->m_abox[0]);
#ifdef _OPENMP
#pragma omp parallel
#endif
	    {
#ifndef _OPENMP
		int tid = 0;
#else
		int tid = omp_get_thread_num();
#pragma omp for
#endif
		for (int i = 0; i < N; ++i) {
		    bxs[tid].minBox(m_ref->m_abox[i]);
		}
	    }
	    minbox = bxs[0];
	    for (int i = 1; i < nthreads; ++i) {
		minbox.minBox(bxs[i]);
	    }
	}
    }
    minbox.coarsen(m_crse_ratio).convert(ixType());
    return minbox;
}

Box
BoxArray::minimalBox (long& npts_avg_box) const
{
    BL_ASSERT(m_simple);
    Box minbox;
    const int N = size();
    long npts_tot = 0;
    if (N > 0)
    {
#ifdef _OPENMP
	bool use_single_thread = omp_in_parallel();
	const int nthreads = use_single_thread ? 1 : omp_get_max_threads();
#else
	bool use_single_thread = true;
	const int nthreads = 1;
#endif
	if (use_single_thread)
	{
	    minbox = m_ref->m_abox[0];
            npts_tot += m_ref->m_abox[0].numPts();
	    for (int i = 1; i < N; ++i) {
		minbox.minBox(m_ref->m_abox[i]);
                npts_tot += m_ref->m_abox[i].numPts();
	    }
	}
	else
	{
	    Vector<Box> bxs(nthreads, m_ref->m_abox[0]);
#ifdef _OPENMP
#pragma omp parallel reduction(+:npts_tot)
#endif
	    {
#ifndef _OPENMP
		int tid = 0;
#else
		int tid = omp_get_thread_num();
#pragma omp for
#endif
		for (int i = 0; i < N; ++i) {
		    bxs[tid].minBox(m_ref->m_abox[i]);
                    long npts = m_ref->m_abox[i].numPts();
                    npts_tot += npts;
		}
	    }
	    minbox = bxs[0];
	    for (int i = 1; i < nthreads; ++i) {
		minbox.minBox(bxs[i]);
	    }
	}
    }
    minbox.coarsen(m_crse_ratio).convert(ixType());
    npts_tot /= AMREX_D_TERM(m_crse_ratio[0],*m_crse_ratio[1],*m_crse_ratio[2]);
    npts_avg_box = npts_tot / N;
    return minbox;
}

bool
BoxArray::intersects (const Box& b, int ng) const
{
    return intersects(b, IntVect(ng));
}

bool
BoxArray::intersects (const Box& b, const IntVect& ng) const
{
    std::vector< std::pair<int,Box> > isects;

    bool first_only = true;
    intersections(b,isects,first_only,ng);

    return (isects.size() > 0) ;
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects,false,IntVect::TheZeroVector());
    return isects;
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx, bool first_only, int ng) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects,first_only,IntVect(ng));
    return isects;
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx, bool first_only, const IntVect& ng) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects,first_only,ng);
    return isects;
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects) const
{
    intersections(bx, isects, false, IntVect::TheZeroVector());
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects,
			 bool                               first_only,
			 int                                ng) const
{
    intersections(bx,isects,first_only,IntVect(ng));
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects,
			 bool                               first_only,
			 const IntVect&                     ng) const
{
  // This is called too many times BL_PROFILE("BoxArray::intersections()");

    BARef::HashType& BoxHashMap = getHashMap();

    isects.resize(0);

    if (!BoxHashMap.empty())
    {
        BL_ASSERT(bx.ixType() == ixType());

	Box gbx = amrex::grow(bx,ng);

	IntVect glo = gbx.smallEnd();
	IntVect ghi = gbx.bigEnd();
	const IntVect& doilo = getDoiLo();
	const IntVect& doihi = getDoiHi();

	gbx.setSmall(glo - doihi).setBig(ghi + doilo);
        gbx.refine(m_crse_ratio).coarsen(m_ref->crsn);
	
        const IntVect& sm = amrex::max(gbx.smallEnd()-1, m_ref->bbox.smallEnd());
        const IntVect& bg = amrex::min(gbx.bigEnd(),     m_ref->bbox.bigEnd());

        Box cbx(sm,bg);
        cbx.normalize();

	if (!cbx.intersects(m_ref->bbox)) return;

	auto TheEnd = BoxHashMap.cend();

        bool super_simple = m_simple && m_crse_ratio==1 && m_typ.cellCentered();
        auto& abox = m_ref->m_abox;

        for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); iv <= End; cbx.next(iv))
        {
            auto it = BoxHashMap.find(iv);

            if (it != TheEnd)
            {
                for (const int index : it->second)
                {
                    const Box& ibox = super_simple ? abox[index] : (*this)[index];
                    const Box& isect = bx & amrex::grow(ibox,ng);

                    if (isect.ok())
                    {
                        isects.push_back(std::pair<int,Box>(index,isect));
			if (first_only) return;
                    }
                }
            }
        }
    }
}

BoxList
BoxArray::complementIn (const Box& bx) const
{
    BoxList bl(bx.ixType());
    complementIn(bl, bx);
    return bl;
}

void
BoxArray::complementIn (BoxList& bl, const Box& bx) const
{
    bl.clear();
    bl.set(bx.ixType());
    bl.push_back(bx);

    if (!empty()) 
    {
	BARef::HashType& BoxHashMap = getHashMap();

	BL_ASSERT(bx.ixType() == ixType());

	Box gbx = bx;

	IntVect glo = gbx.smallEnd();
	IntVect ghi = gbx.bigEnd();
	const IntVect& doilo = getDoiLo();
	const IntVect& doihi = getDoiHi();

	gbx.setSmall(glo - doihi).setBig(ghi + doilo);
        gbx.refine(m_crse_ratio).coarsen(m_ref->crsn);
	
        const IntVect& sm = amrex::max(gbx.smallEnd()-1, m_ref->bbox.smallEnd());
        const IntVect& bg = amrex::min(gbx.bigEnd(),     m_ref->bbox.bigEnd());

        Box cbx(sm,bg);
        cbx.normalize();

	if (!cbx.intersects(m_ref->bbox)) return;

	auto TheEnd = BoxHashMap.cend();

        BoxList newbl(bl.ixType());
        newbl.reserve(bl.capacity());
        BoxList newdiff(bl.ixType());

        bool super_simple = m_simple && m_crse_ratio==1 && m_typ.cellCentered();
        auto& abox = m_ref->m_abox;

	for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); 
	     iv <= End && bl.isNotEmpty(); 
	     cbx.next(iv))
        {
            auto it = BoxHashMap.find(iv);

            if (it != TheEnd)
            {
                for (const int index : it->second)
                {
                    const Box& isect = (super_simple)
                        ? (bx & abox[index])
                        : (bx & (*this)[index]);

                    if (isect.ok())
                    {
                        newbl.clear();
                        for (const Box& b : bl) {
                            amrex::boxDiff(newdiff, b, isect);
                            newbl.join(newdiff);
                        }
                        bl.swap(newbl);
                    }
                }
            }
        }
    }
}

void
BoxArray::clear_hash_bin () const
{
    if (!m_ref->hash.empty())
    {
#ifdef BL_MEM_PROFILING
	m_ref->updateMemoryUsage_hash(-1);
#endif
        m_ref->hash.clear();
        m_ref->has_hashmap = false;
    }
}

//
// Currently this assumes your Boxes are cell-centered.
//
void
BoxArray::removeOverlap (bool simplify)
{
    if (! ixType().cellCentered()) {
        amrex::Abort("BoxArray::removeOverlap() supports cell-centered only");
    }

    if (m_crse_ratio != 1) {
        amrex::Abort("BoxArray::removeOverlap() must have m_crse_ratio == 1");
    }

    uniqify();

    BARef::HashType& BoxHashMap = m_ref->hash;

    const Box EmptyBox;

    std::vector< std::pair<int,Box> > isects;
    //
    // Note that "size()" can increase in this loop!!!
    //
#ifdef BL_MEM_PROFILING
    m_ref->updateMemoryUsage_box(-1);
    m_ref->updateMemoryUsage_hash(-1);
    long total_hash_bytes_save = m_ref->total_hash_bytes;
#endif

    BoxList bl_diff;

    for (int i = 0; i < size(); i++)
    {
        if (m_ref->m_abox[i].ok())
        {
            intersections(m_ref->m_abox[i],isects);

            for (int j = 0, N = isects.size(); j < N; j++)
            {
                if (isects[j].first == i) continue;

                Box& bx = m_ref->m_abox[isects[j].first];

                amrex::boxDiff(bl_diff, bx, isects[j].second);

                bx = EmptyBox;

                for (const Box& b : bl_diff)
                {
                    m_ref->m_abox.push_back(b);
                    BoxHashMap[amrex::coarsen(b.smallEnd(),m_ref->crsn)].push_back(size()-1);
                }
            }
        }
    }
#ifdef BL_MEM_PROFILING
    m_ref->updateMemoryUsage_box(1);
#endif
    //
    // We now have "holes" in our BoxArray. Make us good.
    //
    BoxList bl(ixType());
    for (const auto& b : m_ref->m_abox) {
        if (b.ok()) {
            bl.push_back(b);
        }
    }
    
    if (simplify) {
        bl.simplify();
    }

    BoxArray nba(std::move(bl));

    *this = nba;

#ifdef BL_MEM_PROFILING
    m_ref->total_hash_bytes = total_hash_bytes_save;
#endif

    BL_ASSERT(isDisjoint());
}

void
BoxArray::type_update ()
{
    if (!empty())
    {
	m_typ = m_ref->m_abox[0].ixType();
        m_transformer->setIxType(m_typ);
	if (!m_typ.cellCentered())
	{
            for (auto& bx : m_ref->m_abox) {
		bx.enclosedCells();
	    }
	}
    }
}

IntVect
BoxArray::getDoiLo () const noexcept
{
    return m_simple ? IntVect::TheZeroVector() : m_transformer->doiLo();
}

IntVect
BoxArray::getDoiHi () const noexcept
{
    return m_simple ?           m_typ.ixType() : m_transformer->doiHi();
}

BARef::HashType&
BoxArray::getHashMap () const
{
    BARef::HashType& BoxHashMap = m_ref->hash;

    if (m_ref->HasHashMap()) return BoxHashMap;

#ifdef _OPENMP
    #pragma omp critical(intersections_lock)
#endif
    {
        if (BoxHashMap.empty() && size() > 0)
        {
            //
            // Calculate the bounding box & maximum extent of the boxes.
            //
	    IntVect maxext = IntVect::TheUnitVector();
            Box boundingbox = m_ref->m_abox[0];

	    const int N = size();
	    for (int i = 0; i < N; ++i)
            {
                const Box& bx = m_ref->m_abox[i];
                maxext = amrex::max(maxext, bx.size());
                boundingbox.minBox(bx);
            }

            for (int i = 0; i < N; i++)
            {
                const IntVect& crsnsmlend 
		    = amrex::coarsen(m_ref->m_abox[i].smallEnd(),maxext);
                BoxHashMap[crsnsmlend].push_back(i);
            }

            m_ref->crsn = maxext;
            m_ref->bbox =boundingbox.coarsen(maxext);
            m_ref->bbox.normalize();

	    m_ref->has_hashmap = true;

#ifdef BL_MEM_PROFILING
	    m_ref->updateMemoryUsage_hash(1);
#endif
        }
    }

    return BoxHashMap;
}

void
BoxArray::uniqify ()
{
    if (m_ref.use_count() == 1) {
        clear_hash_bin();
    } else {
	auto p = std::make_shared<BARef>(*m_ref);
	std::swap(m_ref,p);
    }
    if (m_crse_ratio != 1) {
        const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++) {
            m_ref->m_abox[i].coarsen(m_crse_ratio);
        }
        m_crse_ratio = IntVect::TheUnitVector();
        m_transformer->setCrseRatio(m_crse_ratio);
    }
}

std::ostream&
operator<< (std::ostream&   os,
            const BoxArray& ba)
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << "(BoxArray maxbox("
       << ba.size()
       << ")\n       m_ref->m_hash_sig("
       << 0
       << ")\n       ";

    for (int i = 0, N = ba.size(); i < N; ++i)
        os << ba[i] << ' ';

    os << ")\n";

    if (os.fail())
        amrex::Error("operator<<(ostream& os,const BoxArray&) failed");

    return os;
}

BoxArray
boxComplement (const Box& b1in,
               const Box& b2)
{
    return BoxArray(amrex::boxDiff(b1in, b2));
}

BoxArray
complementIn (const Box&      b,
              const BoxArray& ba)
{
    return BoxArray { ba.complementIn(b) };
}

BoxArray
intersect (const BoxArray& ba,
           const Box&      b,
           int   ng)
{
    std::vector< std::pair<int,Box> > isects;

    ba.intersections(b,isects,false,IntVect(ng));

    const int N = isects.size();

    BoxArray r(N);

    if (N > 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
        {
            r.set(i, isects[i].second);
        }
    }

    return r;
}

BoxArray
intersect (const BoxArray& ba,
           const Box&      b,
           const IntVect&  ng)
{
    std::vector< std::pair<int,Box> > isects;

    ba.intersections(b,isects,false,ng);

    const int N = isects.size();

    BoxArray r(N);

    if (N > 0) {
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
        {
            r.set(i, isects[i].second);
        }
    }

    return r;
}

BoxArray
intersect (const BoxArray& lhs, const BoxArray& rhs)
{
    if (lhs.size() == 0 || rhs.size() == 0) return BoxArray();
    BoxList bl(lhs[0].ixType());
    for (int i = 0, Nl = lhs.size(); i < Nl; ++i)
    {
        const BoxArray& ba = amrex::intersect(rhs, lhs[i]);
        bl.join(ba.boxList());
    }
    return BoxArray(bl);
}

BoxList
intersect (const BoxArray& ba, const BoxList& bl)
{
    BoxList newbl(bl.ixType());
    for (const Box& bx : bl)
    {
        const BoxArray& newba = amrex::intersect(ba, bx);
        newbl.join(newba.boxList());
    }
    return newbl;
}

BoxArray
convert (const BoxArray& ba, IndexType typ)
{
    BoxArray ba2 = ba;
    ba2.convert(typ);
    return ba2;
}

BoxArray
convert (const BoxArray& ba, const IntVect& typ)
{
    BoxArray ba2 = ba;
    ba2.convert(IndexType(typ));
    return ba2;
}

BoxArray
coarsen (const BoxArray& ba, int ratio)
{
    BoxArray ba2 = ba;
    ba2.coarsen(ratio);
    return ba2;
}

BoxArray
coarsen (const BoxArray& ba, const IntVect& ratio)
{
    BoxArray ba2 = ba;
    ba2.coarsen(ratio);
    return ba2;
}

BoxList
GetBndryCells (const BoxArray& ba,
               int             ngrow)
{
    BL_ASSERT(ba.ok());
    BL_ASSERT(ba.size() > 0);
    //
    // First get list of all ghost cells.
    //
    const IndexType btype = ba.ixType();

    BoxList bcells = ba.boxList();
    bcells.simplify();

    BoxArray tba(bcells);

    BoxList gcells(btype);
    BoxList bl_diff(btype);
    for (int i = 0, N = tba.size(); i < N; ++i)
    {
	const Box& bx = tba[i];
        amrex::boxDiff(bl_diff, amrex::grow(bx,ngrow), bx);
        gcells.join(bl_diff);
    }
    //
    // Now strip out intersections with original BoxArray.
    //
    std::vector< std::pair<int,Box> > isects;

    bcells.clear();
    BoxList pieces(btype);
    BoxList bl_tmp(btype);

    for (const Box& gbx : gcells)
    {
        tba.intersections(gbx, isects);
        if (isects.empty())
        {
            bcells.push_back(gbx);
        }
        else
        {
            pieces.clear();
            for (const auto& isec : isects) {
                pieces.push_back(isec.second);
            }
            bl_tmp.complementIn(gbx,pieces);
            bcells.join(bl_tmp);
        }
    }

    gcells = amrex::removeOverlap(bcells);
    gcells.simplify();

    return gcells;
}


void
readBoxArray (BoxArray&     ba,
              std::istream& is,
              bool          bReadSpecial)
{
    if (bReadSpecial == false)
    {
        ba.readFrom(is);
    }
    else
    {
        BL_ASSERT(ba.size() == 0);
        int maxbox;
        unsigned long in_hash; // will be ignored
        is.ignore(bl_ignore_max, '(') >> maxbox >> in_hash;
        ba.resize(maxbox);
        for (int i = 0; i < maxbox; i++)
        {
            Box b;
            is >> b;
            ba.set(i, b);
        }
        is.ignore(bl_ignore_max, ')');

        if (is.fail()) {
            amrex::Error("readBoxArray(BoxArray&,istream&,int) failed");
        }
    }
}

bool match (const BoxArray& x, const BoxArray& y)
{
    if (x == y) {
	return true;
    } else {
	bool m = (x.size() == y.size()) && (x.ixType() == y.ixType());
	for (int i = 0, N = x.size(); i < N && m; ++i) {
	    m = x[i] == y[i];
	}
	return m;
    }
}

std::ostream&
operator<< (std::ostream& os, const BoxArray::RefID& id)
{
    os << id.data;
    return os;
}

}
