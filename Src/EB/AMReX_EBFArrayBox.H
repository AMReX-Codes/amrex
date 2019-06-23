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
    EBFArrayBox (Arena* ar);
    EBFArrayBox (const EBCellFlagFab& ebcellflag, const Box& box, int ncomps, Arena* ar);
    EBFArrayBox (EBFArrayBox const& rhs, MakeType make_type, int scomp, int ncomp);

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    ~EBFArrayBox ();

    const EBCellFlagFab& getEBCellFlagFab () const { return *m_ebcellflag; }

private:
    const EBCellFlagFab* m_ebcellflag;
};

const EBCellFlagFab& getEBCellFlagFab (const FArrayBox& fab);

}

#endif
