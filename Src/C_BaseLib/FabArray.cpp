#include <FabArray.H>

FabArrayBase::FabArrayBase()
{
}

FabArrayBase::FabArrayBase(const BoxArray& bxs, int nvar, int ngrow)
    :
    boxarray(bxs),
    distributionMap(boxarray, ParallelDescriptor::NProcsCFD()),
    n_grow(ngrow),
    n_comp(nvar)
{
}

FabArrayBase::FabArrayBase(const BoxArray& bxs, int nvar, int ngrow, const DistributionMapping& map)
    :
    boxarray(bxs),
    distributionMap(map),
    n_grow(ngrow),
    n_comp(nvar)
{
}

FabArrayBase::~FabArrayBase()
{
}

int
FabArrayBase::nGrow () const
{
    return n_grow;
}

const BoxArray&
FabArrayBase::boxArray () const
{
    return boxarray;
}

const Box&
FabArrayBase::box (int K) const
{
    return boxarray[K];
}

Box
FabArrayBase::fabbox (int K) const
{
    //
    // Do not use fabparray[K] because it may not be valid in parallel.
    //
    return BoxLib::grow(boxarray[K], n_grow);
}

size_t
FabArrayBase::size () const
{
    return boxarray.size();
}

int
FabArrayBase::nComp () const
{
    return n_comp;
}

const DistributionMapping&
FabArrayBase::DistributionMap () const
{
    return distributionMap;
}

MFIter::MFIter (const FabArrayBase& fabarray)
    :
    fabArray(fabarray),
    currentIndex(0)
{
#ifdef BL_USE_MPI
    //
    // Increment the currentIndex to start at the first valid index
    // for this ParallelDescriptor::MyProc.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    while (fabArray.DistributionMap()[currentIndex] != MyProc)
    {
        ++currentIndex;
    }
#endif
}

MFIter&
MFIter::operator++ ()
{
#ifdef BL_USE_MPI
    const int MyProc = ParallelDescriptor::MyProc();
    //
    // Go to the next index on this processor.
    //
    do
    {
        ++currentIndex;
    }
    while (fabArray.DistributionMap()[currentIndex] != MyProc);
#else
    ++currentIndex;
#endif

    return *this;
}

MFIter
MFIter::operator++ (int)
{
    MFIter res(*this);
    this->operator++();
    return res;
}

int
MFIter::index () const
{
    return currentIndex;
}

const Box&
MFIter::validbox () const
{
    return fabArray.box(currentIndex);
}

Box
MFIter::fabbox () const
{
    return fabArray.fabbox(currentIndex);
}

bool
MFIter::isValid ()
{
    BL_ASSERT(currentIndex >= 0);

#if defined(BL_USE_MPI) && defined(HG_DEBUG)
    extern bool HG_is_debugging;
    if (HG_is_debugging)
    {
	if (currentIndex < fabArray.size()) return true;
	BL_MPI_REQUIRE( MPI_Barrier(ParallelDescriptor::Communicator()) );
	return false;
    }
#endif
    return currentIndex < fabArray.size();
}
