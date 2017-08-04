#ifndef AMREX_EBFARRAYBOX_H_
#define AMREX_EBFARRAYBOX_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBISBox.H>

namespace amrex {

class EBFlagFab;

class EBFArrayBox
    : public FArrayBox
{
public:
    EBFArrayBox ();
    EBFArrayBox (const EBISBox& ebisBox, const EBFlagFab& ebflag, const Box& box, int ncomps);

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    virtual ~EBFArrayBox ();

private:
    EBISBox m_ebisbox;
};

}

#endif
