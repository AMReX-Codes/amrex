#ifndef AMREX_EBFARRAYBOX_H_
#define AMREX_EBFARRAYBOX_H_
#include <AMReX_Config.H>

#include <AMReX_FArrayBox.H>

namespace amrex {

class EBCellFlagFab;
class EBFArrayBoxFactory;

class EBFArrayBox
    : public FArrayBox
{
public:
    EBFArrayBox ();
    explicit EBFArrayBox (Arena* ar);
    EBFArrayBox (const EBCellFlagFab& ebflag, const Box& bx, int ncomps, Arena* ar,
                 const EBFArrayBoxFactory* factory=nullptr, int box_index=-1);
    EBFArrayBox (EBFArrayBox const& rhs, MakeType make_type, int scomp, int ncomp);

    EBFArrayBox (EBFArrayBox&& rhs) noexcept = default;
    EBFArrayBox (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (const EBFArrayBox&) = delete;
    EBFArrayBox& operator= (EBFArrayBox&&) = delete;

    ~EBFArrayBox () override;

    //! Get EBCellFlag Fab
    [[nodiscard]] const EBCellFlagFab& getEBCellFlagFab () const { return *m_ebcellflag; }

    //! Get a pointer to levelset data if available.  The return value could
    //! be nullptr if not available.
    [[nodiscard]] const FArrayBox* getLevelSetData () const;

    //! Get a pointer to volume fraction data if available.  The return
    //! value could be nullptr if not available.
    [[nodiscard]] const FArrayBox* getVolFracData () const;

    //! Get a pointer to volume centroid data if available.  The return
    //! value could be nullptr if not available.
    [[nodiscard]] const FArrayBox* getCentroidData () const;

    //! Get a pointer to boundary centroid data if available.  The return
    //! value could be nullptr if not available.
    [[nodiscard]] const FArrayBox* getBndryCentData () const;

    //! Get a pointer to boundary normal data if available.  The return
    //! value could be nullptr if not available.
    [[nodiscard]] const FArrayBox* getBndryNormalData () const;

    //! Get a pointer to boundary area data if available.  The return
    //! value could be nullptr if not available.
    [[nodiscard]] const FArrayBox* getBndryAreaData () const;

    //! Get pointers to area fraction data if available.  The return value
    //! could be nullptr if not available.
    [[nodiscard]] Array<const FArrayBox*, AMREX_SPACEDIM> getAreaFracData () const;

    //! Get pointers to face centroid data if available.  The return value
    //! could be nullptr if not available.
    [[nodiscard]] Array<const FArrayBox*, AMREX_SPACEDIM> getFaceCentData () const;

    //! Get pointers to edge centroid data if available.  The return value
    //! could be nullptr if not available.
    [[nodiscard]] Array<const FArrayBox*, AMREX_SPACEDIM> getEdgeCentData () const;

private:
    const EBCellFlagFab* m_ebcellflag = nullptr;
    const EBFArrayBoxFactory* m_factory = nullptr;
    int m_box_index = -1;
};

const EBCellFlagFab& getEBCellFlagFab (const FArrayBox& fab);

}

#endif
