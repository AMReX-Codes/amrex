
#include <AMReX_IndexType.H>

#include <iostream>
#include <iomanip>

namespace amrex::detail {

std::ostream&
index_type_write (std::ostream& os, const unsigned int& iv, int dim)
{
    os << '(' << ((iv & 1u != 0) ? 'N' : 'C');
    for (int i=1; i<dim; ++i) {
        os << ',' << ((iv & (1u<<i) != 0) ? 'N' : 'C');
    }
    os << ')' << std::flush;

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
index_type_read (std::istream& is, unsigned int& iv, int dim)
{
    char t = '0';
    is.ignore(BL_IGNORE_MAX, '(') >> t;
    BL_ASSERT(t == 'C' || t == 'N');
    t == 'N' ? (iv |= 1u) : (iv &= ~1u);
    for (int i=1; i<dim; ++i) {
        is.ignore(BL_IGNORE_MAX, ',') >> t;
        BL_ASSERT(t == 'C' || t == 'N');
        t == 'N' ? (iv |= (1u << i)) : (iv &= ~(1u << i));
    }
    is.ignore(BL_IGNORE_MAX, ')');

    if (is.fail()) {
        amrex::Error("operator>>(ostream&,IndexType&) failed");
    }

    return is;
}

}
