
#include <AMReX_IndexType.H>

#include <iostream>
#include <iomanip>

namespace amrex {

std::ostream&
operator<< (std::ostream&    os,
            const IndexType& it)
{
    os << '('
       << AMREX_D_TERM( (it.test(0)?'N':'C'),
                  << ',' << (it.test(1)?'N':'C'),
                  << ',' << (it.test(2)?'N':'C')) << ')' << std::flush;

    if (os.fail()) {
        amrex::Error("operator<<(ostream&,IndexType&) failed");
    }

    return os;
}

//
// Copied from <Utility.H>
//
#define BL_IGNORE_MAX 100000

std::istream&
operator>> (std::istream& is,
            IndexType&    it)
{
    char AMREX_D_DECL(t0,t1,t2);

    AMREX_D_EXPR( is.ignore(BL_IGNORE_MAX, '(') >> t0,
            is.ignore(BL_IGNORE_MAX, ',') >> t1,
            is.ignore(BL_IGNORE_MAX, ',') >> t2);
    is.ignore(BL_IGNORE_MAX, ')');
    AMREX_D_TERM(
        BL_ASSERT(t0 == 'C' || t0 == 'N'); t0=='N'?it.set(0):it.unset(0); ,
        BL_ASSERT(t1 == 'C' || t1 == 'N'); t1=='N'?it.set(1):it.unset(1); ,
        BL_ASSERT(t2 == 'C' || t2 == 'N'); t2=='N'?it.set(2):it.unset(2));

    if (is.fail()) {
        amrex::Error("operator>>(ostream&,IndexType&) failed");
    }

    return is;
}

// IndexTypeSet::IndexTypeSet (unsigned long a_iTypeSet) : m_EncodedLong(a_iTypeSet) {};

// IndexTypeSet::IndexTypeSet (Vector<IndexType> a_IxTypes)
// {
//     define(a_IxTypes);
// }

// IndexTypeSet::IndexTypeSet (Vector<IntVect> a_IntVects)
// {
//     define(ConvertToIxTypes(a_IntVects));
// }

// void
// IndexTypeSet::define(Vector<IndexType> a_IxTypes)
// {
//     m_n = a_IxTypes.size();
//     m_EncodedLong = ConvertToEncodedLong(a_IxTypes);
// }

// Vector<IndexType>
// IndexTypeSet::IndexTypeList ()
// {
//     return ConvertToIxTypes(m_EncodedLong);
// }

// Vector<IntVect>
// IndexTypeSet::IntVectList ()
// {
//     return ConvertToIntVects(m_EncodedLong);
// }

constexpr std::size_t
IndexTypeSet::GetSize (unsigned long a_EncodedLong)
{
    return 1UL + (a_EncodedLong & (IndexTypeSet::MAX_SIZE - 1UL));
}

template <std::size_t N>
constexpr Array<IndexType,N>
IndexTypeSet::ConvertToIxTypes (unsigned long a_EncodedLong)
{   
    BL_ASSERT(N == IndexTypeSet::GetSize(a_EncodedLong));
    Array<IndexType,N> o_IxTypes;
    auto tmp = a_EncodedLong >> IndexTypeSet::NUM_BITS_FOR_SIZE;
    for (int i = N-1; i >= 0; ++i) {
        IndexType ixtype = IndexType();
        ixtype.itype = tmp & (1 << AMREX_SPACEDIM);
        o_IxTypes[i] = ixtype;
        tmp = tmp >> AMREX_SPACEDIM;
    }
    return o_IxTypes;
}

template <std::size_t N>
constexpr Array<IntVect,N>
IndexTypeSet::ConvertToIntVects (unsigned long a_EncodedLong)
{
    return ConvertToIntVects(ConvertToIxTypes(a_EncodedLong));
}

template <std::size_t N>
constexpr unsigned long
IndexTypeSet::ConvertToEncodedLong (Array<IndexType,N> a_IxTypes)
{
    BL_ASSERT(N <= IndexTypeSet::MAX_SIZE);
    unsigned long o_EncodedLong = 0;
    for (IndexType ixtype : a_IxTypes) {
        o_EncodedLong = o_EncodedLong << AMREX_SPACEDIM;
        o_EncodedLong += ixtype.itype;
    }
    o_EncodedLong = o_EncodedLong << IndexTypeSet::NUM_BITS_FOR_SIZE;
    o_EncodedLong += N;
    return o_EncodedLong;
}

template <std::size_t N>
constexpr unsigned long
IndexTypeSet::ConvertToEncodedLong (Array<IntVect,N> a_IntVects)
{
    return ConvertToEncodedLong(ConvertToIxTypes(a_IntVects));
}

template <std::size_t N>
constexpr Array<IndexType,N>
IndexTypeSet::ConvertToIxTypes (Array<IntVect,N> a_IntVects)
{
    Array<IndexType,N> o_IxTypes;
    for (int i = 0; i < N; ++i) {
        o_IxTypes[i] = IndexType(a_IntVects[i]);
    }
    return o_IxTypes;
}

template <std::size_t N>
constexpr Array<IntVect,N>
IndexTypeSet::ConvertToIntVects (Array<IndexType,N> a_IxTypes)
{
    Array<IntVect,N> o_IntVects;
    for (int i = 0; i < N; ++i) {
        o_IntVects[i] = a_IxTypes[i].toIntVect();
    }
    return o_IntVects;
}

}
