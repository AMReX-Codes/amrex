
#include <AMReX_BLassert.H>
#include <AMReX_BoxArray.H>
#include <AMReX_ParallelDescriptor.H>
#include <AMReX_Utility.H>

#ifdef BL_MEM_PROFILING
#include <AMReX_MemProfiler.H>
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
{ 
    define(bl); 
}

BARef::BARef (std::istream& is)
{ 
    define(is); 
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

ptrdiff_t 
BARef::getRefID () const
{
    static BARef ref0;
    return &(*this) - &(ref0);
}

void
BARef::define (std::istream& is)
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    BL_ASSERT(m_abox.size() == 0);
    int           maxbox;
    unsigned long hash;
    is.ignore(bl_ignore_max, '(') >> maxbox >> hash;
    resize(maxbox);
    for (Array<Box>::iterator it = m_abox.begin(), End = m_abox.end(); it != End; ++it)
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
    BL_ASSERT(m_abox.size() == 0);
    const int N = bl.size();
    resize(N);
    int count = 0;
    for (BoxList::const_iterator bli = bl.begin(), End = bl.end(); bli != End; ++bli)
    {
        m_abox[count++] = *bli;
    }
}

void 
BARef::resize (long n) {
#ifdef BL_MEM_PROFILING
    updateMemoryUsage_box(-1);
    updateMemoryUsage_hash(-1);
#endif
    m_abox.resize(n);
    hash.clear();
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
}

void
BoxArray::Initialize ()
{
    if (!initialized) {
	initialized = true;
	BARef::Initialize();
    }
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
    m_crse_ratio(trans.crseRatio()),
    m_simple(trans.simple()),
    m_ref(rhs.m_ref)
{}

BoxArray::BoxArray (const BoxArray& rhs)
    :
    m_transformer(rhs.m_transformer->clone()),
    m_typ(rhs.m_typ),
    m_crse_ratio(rhs.m_crse_ratio),
    m_simple(rhs.m_simple),
    m_ref(rhs.m_ref)
{}

BoxArray::BoxArray(BoxArray&& rhs) noexcept
    :
    m_transformer(std::move(rhs.m_transformer)),
    m_typ(rhs.m_typ),
    m_crse_ratio(rhs.m_crse_ratio),
    m_simple(rhs.m_simple),
    m_ref(std::move(rhs.m_ref))
{}

BoxArray::~BoxArray ()
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
BoxArray::clear ()
{
    m_transformer.reset(new DefaultBATransformer());
    m_ref.reset(new BARef());
}

void
BoxArray::resize (long len)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    m_ref->resize(len);
}

long
BoxArray::size () const
{
    return m_ref->m_abox.size();
}

long
BoxArray::capacity () const
{
    return m_ref->m_abox.capacity();
}

bool
BoxArray::empty () const
{
    return m_ref->m_abox.empty();
}

long
BoxArray::numPts () const
{
    long result = 0;
    const int N = size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
    for (int i = 0; i < N; ++i)
    {
        result += get(i).numPts();
    }
    return result;
}

double
BoxArray::d_numPts () const
{
    double result = 0;
    const int N = size();
#ifdef _OPENMP
#pragma omp parallel for reduction(+:result)
#endif
    for (int i = 0; i < N; ++i)
    {
        result += get(i).d_numPts();
    }
    return result;
}

void
BoxArray::readFrom (std::istream& is)
{
    BL_ASSERT(size() == 0);
    clear();
    m_simple = true;
    m_crse_ratio = IntVect::TheUnitVector();
    m_ref->define(is);
    type_update();
}

std::ostream&
BoxArray::writeOn (std::ostream& os) const
{
    //
    // TODO -- completely remove the fiction of a hash value.
    //
    os << '(' << size() << ' ' << 0 << '\n';

    const int N = size();
    for (int i = 0; i < N; ++i)
        os << get(i) << '\n';

    os << ')';

    if (os.fail())
        amrex::Error("BoxArray::writeOn(ostream&) failed");

    return os;
}

bool
BoxArray::operator== (const BoxArray& rhs) const
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
BoxArray::operator!= (const BoxArray& rhs) const
{
    return !operator==(rhs);
}

bool
BoxArray::CellEqual (const BoxArray& rhs) const
{
    return m_crse_ratio == rhs.m_crse_ratio
        && (m_ref == rhs.m_ref || m_ref->m_abox == rhs.m_ref->m_abox);
}

BoxArray&
BoxArray::maxSize (int block_size)
{
    return maxSize(IntVect(D_DECL(block_size,block_size,block_size)));
}

BoxArray&
BoxArray::maxSize (const IntVect& block_size)
{
    BoxList blst(*this);
    blst.maxSize(block_size);
    const int N = blst.size();
    if (size() != N) { // If size doesn't change, do nothing.
	m_ref = std::make_shared<BARef>(blst);
    }
    return *this;
}

BoxArray&
BoxArray::refine (int refinement_ratio)
{
    return refine(IntVect(D_DECL(refinement_ratio,refinement_ratio,refinement_ratio)));
}

BoxArray&
BoxArray::refine (const IntVect& iv)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
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

BoxArray&
BoxArray::coarsen (int refinement_ratio)
{
    return coarsen(IntVect(D_DECL(refinement_ratio,refinement_ratio,refinement_ratio)));
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
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(n).coarsen(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int n)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(n);
    return *this;
}

BoxArray&
BoxArray::grow (const IntVect& iv)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(iv);
    return *this;
}

BoxArray&
BoxArray::grow (int dir,
                int n_cell)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].grow(dir, n_cell);
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
BoxArray::convert (Box (*fp)(const Box&))
{
    BL_ASSERT(!(fp == 0));

    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = size();
    for (int i = 0; i < N; ++i)
	set(i,fp(get(i)));
    return *this;
}

BoxArray&
BoxArray::shift (int dir,
                 int nzones)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shift(dir, nzones);
    return *this;
}

BoxArray&
BoxArray::shift (const IntVect& iv)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shift(iv);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (int dir,
                     int num_halfs)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shiftHalf(dir, num_halfs);
    return *this;
}

BoxArray&
BoxArray::shiftHalf (const IntVect& iv)
{
    if (m_ref.use_count()==1) {
	clear_hash_bin();
    } else {
        uniqify();
    }
    const int N = m_ref->m_abox.size();
#ifdef _OPENMP
#pragma omp parallel for
#endif
    for (int i = 0; i < N; i++)
        m_ref->m_abox[i].shiftHalf(iv);
    return *this;
}

void
BoxArray::set (int        i,
               const Box& ibox)
{
    if (i == 0) {
	m_transformer->setIxType(ibox.ixType());
        m_typ = ibox.ixType();

	if (m_ref.use_count()==1) {
	    clear_hash_bin();
	} else {
	    uniqify();
	}
    }

    m_ref->m_abox[i] = amrex::enclosedCells(ibox);
}

bool
BoxArray::ok () const
{
    bool isok = true;

    if (size() > 0)
    {
        if (size() == 1) isok = get(0).ok();

        for (int i = 1, N = size(); i < N && isok; i++)
        {
	    isok = get(i).ok();
        }
    }

    return isok;
}

bool
BoxArray::isDisjoint () const
{
    BL_ASSERT(ixType().cellCentered());

    std::vector< std::pair<int,Box> > isects;

    const int N = size();
    for (int i = 0; i < N; ++i)
    {
	intersections(get(i),isects);
        if ( isects.size() > 1 )
            return false;
    }

    return true;
}

BoxList
BoxArray::boxList () const
{
    BoxList newb;
    const int N = size();
    if ( N > 0 ) {
	newb.set(ixType());
	for (int i = 0; i < N; ++i)
	    newb.push_back(get(i));
    }
    return newb;
}

bool
BoxArray::contains (const IntVect& iv) const
{
    if (size() > 0)
    {
	return intersects(Box(iv,iv,ixType()));
    } else {
	return false;
    }
}

bool
BoxArray::contains (const Box& b, bool assume_disjoint_ba) const
{
    BL_ASSERT(!assume_disjoint_ba || isDisjoint());

    bool result = false;

    if (size() > 0)
    {
        BL_ASSERT(ixType() == b.ixType());

        std::vector< std::pair<int,Box> > isects;

        intersections(b,isects);

        if (isects.size() > 0)
        {
	    if (assume_disjoint_ba) {
		long nbx = b.numPts(), nisects = 0L;
		for (int i = 0, N = isects.size(); i < N; i++) {
		    nisects += isects[i].second.numPts();
		}
		result = nbx == nisects;
	    } else {
		BoxList bl(b.ixType());
		for (int i = 0, N = isects.size(); i < N; i++) {
		    bl.push_back(isects[i].second);
		}
		result = bl.contains(b);
	    }
        }
    }

    return result;
}

bool
BoxArray::contains (const BoxArray& bl, bool assume_disjoint_ba) const
{
    if (size() == 0) return false;

    if (!minimalBox().contains(bl.minimalBox())) return false;

    for (int i = 0, N = bl.size(); i < N; ++i) {
        if (!contains(bl[i],assume_disjoint_ba)) {
            return false;
	}
    }

    return true;
}

Box
BoxArray::minimalBox () const
{
    Box minbox;
    const int N = size();
    if (N > 0)
    {
        minbox = m_ref->m_abox[0];
	for (int i = 1; i < N; ++i)
            minbox.minBox(m_ref->m_abox[i]);
    }
    minbox.convert(ixType());
    return minbox;
}

bool
BoxArray::intersects (const Box& b, int ng) const
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
    intersections(bx,isects,false,0);
    return isects;
}

std::vector< std::pair<int,Box> >
BoxArray::intersections (const Box& bx, bool first_only, int ng) const
{
    std::vector< std::pair<int,Box> > isects;
    intersections(bx,isects,first_only,ng);
    return isects;
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects) const
{
    intersections(bx, isects, false, 0);
}

void
BoxArray::intersections (const Box&                         bx,
                         std::vector< std::pair<int,Box> >& isects,
			 bool                               first_only,
			 int                                ng) const
{
    // called too many times  BL_PROFILE("BoxArray::intersections()");

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
        gbx.coarsen(m_ref->crsn);
	
        const IntVect& sm = amrex::max(gbx.smallEnd()-1, m_ref->bbox.smallEnd());
        const IntVect& bg = amrex::min(gbx.bigEnd(),     m_ref->bbox.bigEnd());

        Box cbx(sm,bg);

	if (!cbx.intersects(m_ref->bbox)) return;

	BARef::HashType::const_iterator TheEnd = BoxHashMap.end();

        for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); iv <= End; cbx.next(iv))
        {
            BARef::HashType::const_iterator it = BoxHashMap.find(iv);

            if (it != TheEnd)
            {
                for (std::vector<int>::const_iterator v_it = it->second.begin(), v_end = it->second.end();
                     v_it != v_end;
                     ++v_it)
                {
                    const int  index = *v_it;
                    const Box& isect = bx & amrex::grow(get(index),ng);

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
BoxArray::complement (const Box& bx) const
{
    BoxList bl(bx);

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
        gbx.coarsen(m_ref->crsn);
	
        const IntVect& sm = amrex::max(gbx.smallEnd()-1, m_ref->bbox.smallEnd());
        const IntVect& bg = amrex::min(gbx.bigEnd(),     m_ref->bbox.bigEnd());

        Box cbx(sm,bg);

	if (!cbx.intersects(m_ref->bbox)) return bl;

	BARef::HashType::const_iterator TheEnd = BoxHashMap.end();

	for (IntVect iv = cbx.smallEnd(), End = cbx.bigEnd(); 
	     iv <= End && bl.isNotEmpty(); 
	     cbx.next(iv))
        {
            BARef::HashType::const_iterator it = BoxHashMap.find(iv);

            if (it != TheEnd)
            {
                for (std::vector<int>::const_iterator v_it = it->second.begin(), v_end = it->second.end(); 
		     v_it != v_end && bl.isNotEmpty(); 
		     ++v_it)
                {
                    const int  index = *v_it;
                    const Box& isect = bx & get(index);

                    if (isect.ok())
                    {
			for (BoxList::iterator bli = bl.begin(); bli != bl.end(); )
			{
			    BoxList diff = amrex::boxDiff(*bli, isect);
			    bl.splice_front(diff);
			    bl.remove(bli++);
			}
                    }
                }
            }
        }
    }

    return bl;
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
    }
}

//
// Currently this assumes your Boxes are cell-centered.
//
void
BoxArray::removeOverlap ()
{
    BL_ASSERT(ixType().cellCentered());

    if (m_ref.use_count() > 1) {
        uniqify();
    }

    BARef::HashType& BoxHashMap = m_ref->hash;

    BoxList bl;

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
    for (int i = 0; i < size(); i++)
    {
        if (m_ref->m_abox[i].ok())
        {
            intersections(m_ref->m_abox[i],isects);

            for (int j = 0, N = isects.size(); j < N; j++)
            {
                if (isects[j].first == i) continue;

                Box& bx = m_ref->m_abox[isects[j].first];

                bl = amrex::boxDiff(bx, isects[j].second);

                bx = EmptyBox;

                for (BoxList::const_iterator it = bl.begin(), End = bl.end(); it != End; ++it)
                {
                    m_ref->m_abox.push_back(*it);

                    BoxHashMap[amrex::coarsen(it->smallEnd(),m_ref->crsn)].push_back(size()-1);
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
    bl.clear();

    const Box& bb = m_ref->bbox;

    BARef::HashType::const_iterator TheEnd = BoxHashMap.end();

    for (IntVect iv = bb.smallEnd(), End = bb.bigEnd(); iv <= End; bb.next(iv))
    {
        BARef::HashType::const_iterator it = BoxHashMap.find(iv);

        if (it != TheEnd)
        {
            for (std::vector<int>::const_iterator v_it = it->second.begin(), v_end = it->second.end();
                 v_it != v_end;
                 ++v_it)
            {
                const int index = *v_it;

                if (m_ref->m_abox[index].ok())
                {
                    bl.push_back(m_ref->m_abox[index]);
                }
            }
        }
    }

    bl.simplify();

    BoxArray nba(bl);

    *this = nba;

#ifdef BL_MEM_PROFILING
    m_ref->total_hash_bytes = total_hash_bytes_save;
#endif

    BL_ASSERT(isDisjoint());
}

#ifdef BL_USE_MPI
void
BoxArray::SendBoxArray(const BoxArray &ba, int whichSidecar)
{
    const int MPI_IntraGroup_Broadcast_Rank = ParallelDescriptor::IOProcessor() ? MPI_ROOT : MPI_PROC_NULL;
    Array<int> ba_serial = amrex::SerializeBoxArray(ba);
    int ba_serial_size = ba_serial.size();
    ParallelDescriptor::Bcast(&ba_serial_size, 1, MPI_IntraGroup_Broadcast_Rank,
                              ParallelDescriptor::CommunicatorInter(whichSidecar));
    ParallelDescriptor::Bcast(ba_serial.dataPtr(), ba_serial_size, MPI_IntraGroup_Broadcast_Rank,
                              ParallelDescriptor::CommunicatorInter(whichSidecar));
}

void
BoxArray::RecvBoxArray(BoxArray &ba, int whichSidecar)
{
    int ba_serial_size;
    ParallelDescriptor::Bcast(&ba_serial_size, 1, 0,
                              ParallelDescriptor::CommunicatorInter(whichSidecar));
    Array<int> ba_serial(ba_serial_size);
    ParallelDescriptor::Bcast(ba_serial.dataPtr(), ba_serial_size, 0,
                              ParallelDescriptor::CommunicatorInter(whichSidecar));
    ba = amrex::UnSerializeBoxArray(ba_serial);
}
#endif

void
BoxArray::type_update ()
{
    if (!empty())
    {
	m_typ = m_ref->m_abox[0].ixType();
        m_transformer->setIxType(m_typ);
	if (!m_typ.cellCentered())
	{
	    for (Array<Box>::iterator it = m_ref->m_abox.begin(), End = m_ref->m_abox.end(); 
		 it != End; ++it)
	    {
		it->enclosedCells();
	    }
	}
    }
}

BARef::HashType&
BoxArray::getHashMap () const
{
    BARef::HashType& BoxHashMap = m_ref->hash;

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

	    const int N = size();

	    for (int i = 0; i < N; ++i)
            {
                const Box& bx = m_ref->m_abox[i];
                maxext = amrex::max(maxext, bx.size());
            }

	    IntVect lob = IntVect::TheMaxVector();
	    IntVect hib = IntVect::TheMinVector();

            for (int i = 0; i < N; i++)
            {
                const IntVect& crsnsmlend 
		    = amrex::coarsen(m_ref->m_abox[i].smallEnd(),maxext);
                BoxHashMap[crsnsmlend].push_back(i);
		lob = amrex::min(lob, crsnsmlend);
		hib = amrex::max(hib, crsnsmlend);
            }

            m_ref->crsn = maxext;
            m_ref->bbox = Box(lob,hib);
	    
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
    if (m_ref.use_count() > 1) {
	auto p = std::make_shared<BARef>(*m_ref);
	std::swap(m_ref,p);
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
    return BoxArray(amrex::complementIn(b,ba.boxList()));
}

BoxArray
intersect (const BoxArray& ba,
		   const Box&      b,
		   int   ng)
{
    std::vector< std::pair<int,Box> > isects;

    ba.intersections(b,isects,false,ng);

    BoxArray r(isects.size());

    for (int i = 0, N = isects.size(); i < N; i++)
    {
        r.set(i, isects[i].second);
    }

    return r;
}

BoxArray
intersect (const BoxArray& lhs,
		   const BoxArray& rhs)
{
    if (lhs.size() == 0 || rhs.size() == 0) return BoxArray();
    BoxList bl(lhs[0].ixType());
    for (int i = 0, Nl = lhs.size(); i < Nl; ++i)
    {
        BoxArray ba  = amrex::intersect(rhs, lhs[i]);
        BoxList  tmp = ba.boxList();
        bl.catenate(tmp);
    }
    return BoxArray(bl);
}

BoxArray
convert (const BoxArray& ba, IndexType typ)
{
    BoxArray ba2 = ba;
    return ba2.convert(typ);
}

BoxArray
convert (const BoxArray& ba, const IntVect& typ)
{
    BoxArray ba2 = ba;
    return ba2.convert(IndexType(typ));
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

    BoxList bcells = ba.boxList(), gcells(btype);

    bcells.simplify();

    BoxArray tba(bcells);

    bcells.clear();

    for (int i = 0, N = tba.size(); i < N; ++i)
    {
	const Box& bx = tba[i];
        bcells = amrex::boxDiff(amrex::grow(bx,ngrow),bx);

	gcells.catenate(bcells);
    }
    //
    // Now strip out intersections with original BoxArray.
    //
    std::vector< std::pair<int,Box> > isects;

    for (BoxList::const_iterator it = gcells.begin(), End = gcells.end(); it != End; ++it)
    {
        tba.intersections(*it,isects);

        if (isects.empty())
        {
            bcells.push_back(*it);
        }
        else
        {
            //
            // Collect all the intersection pieces.
            //
            BoxList pieces(btype), leftover(btype);
            for (int i = 0, N = isects.size(); i < N; i++)
                pieces.push_back(isects[i].second);
            leftover = amrex::complementIn(*it,pieces);
            bcells.catenate(leftover);
        }
    }

    //
    // Now strip out overlaps.
    //
    gcells.clear();

    gcells = amrex::removeOverlap(bcells);

    bcells.clear();

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

        if (is.fail())
            amrex::Error("readBoxArray(BoxArray&,istream&,int) failed");
    }
}


Array<int> SerializeBoxArray(const BoxArray &ba)
{
  int nIntsInBox(3 * BL_SPACEDIM);
  Array<int> retArray(ba.size() * nIntsInBox, -1);
  for(int i(0); i < ba.size(); ++i) {
    Array<int> aiBox(amrex::SerializeBox(ba[i]));
    BL_ASSERT(aiBox.size() == nIntsInBox);
    for(int j(0); j < aiBox.size(); ++j) {
      retArray[i * aiBox.size() + j] = aiBox[j];
    }
  }
  return retArray;
}


BoxArray UnSerializeBoxArray(const Array<int> &serarray)
{
  int nIntsInBox(3 * BL_SPACEDIM);
  int nBoxes(serarray.size() / nIntsInBox);
  BoxArray ba(nBoxes);
  for(int i(0); i < nBoxes; ++i) {
    Array<int> aiBox(nIntsInBox);
    for(int j(0); j < nIntsInBox; ++j) {
      aiBox[j] = serarray[i * nIntsInBox + j];
    }
    ba.set(i, amrex::UnSerializeBox(aiBox));
  }
  return ba;
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

}
