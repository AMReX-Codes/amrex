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

    const EBISBox& getEBISBox () const { return m_ebisbox; }
    const EBFlagFab& getEBFlagFab () const { return *m_ebflag; }

private:
    EBISBox m_ebisbox;             // We have to make a copy because EBISL returns a copy
    const EBFlagFab* m_ebflag;
};

const EBFlagFab& getEBFlagFab (const FArrayBox& fab);

const EBISBox& getEBISBox (const FArrayBox& fab);

}

#endif
