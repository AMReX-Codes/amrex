#ifndef AMREX_EBFARRAYBOX_H_
#define AMREX_EBFARRAYBOX_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

class EBFArrayBox
    : public FArrayBox
{
    EBFArrayBox ();

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    virtual ~EBFArrayBox ();
};

}

#endif
