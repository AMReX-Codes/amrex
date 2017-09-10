#ifndef AMREX_EBFARRAYBOX_H_
#define AMREX_EBFARRAYBOX_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_EBISBox.H>

namespace amrex {

class EBCellFlagFab;

class EBFArrayBox
    : public FArrayBox
{
public:
    EBFArrayBox ();
    EBFArrayBox (const EBISBox& ebisBox, const EBCellFlagFab& ebcellflag, const Box& box, int ncomps);

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    virtual ~EBFArrayBox ();

    const EBISBox& getEBISBox () const { return m_ebisbox; }
    const EBCellFlagFab& getEBCellFlagFab () const { return *m_ebcellflag; }

private:
    EBISBox m_ebisbox;             // We have to make a copy because EBISL returns a copy
    const EBCellFlagFab* m_ebcellflag;
};

const EBCellFlagFab& getEBCellFlagFab (const FArrayBox& fab);

const EBISBox& getEBISBox (const FArrayBox& fab);

}

#endif
