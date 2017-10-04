#ifndef AMREX_EBFARRAYBOX_H_
#define AMREX_EBFARRAYBOX_H_

#include <AMReX_FArrayBox.H>

namespace amrex {

class EBCellFlagFab;

class EBFArrayBox
    : public FArrayBox
{
public:
    EBFArrayBox ();
    EBFArrayBox (const EBCellFlagFab& ebcellflag, const Box& box, int ncomps);

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    virtual ~EBFArrayBox ();

    const EBCellFlagFab& getEBCellFlagFab () const { return *m_ebcellflag; }

private:
    const EBCellFlagFab* m_ebcellflag;
};

const EBCellFlagFab& getEBCellFlagFab (const FArrayBox& fab);

}

#endif
