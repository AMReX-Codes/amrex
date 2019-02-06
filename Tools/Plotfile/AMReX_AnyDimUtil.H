#ifndef AMREX_ANY_DIM_UTIL_H_
#define AMREX_ANY_DIM_UTIL_H_

#include <AMReX_Box.H>
#include <iostream>

namespace amrex {

class BoxND
{
public:
    friend std::ostream& operator<< (std::ostream& os, const BoxND& b);
    BoxND (Box const& b, int dim) : m_box(b), m_dim(dim) {}
private:
    Box m_box;
    int m_dim;
};

}

#endif
