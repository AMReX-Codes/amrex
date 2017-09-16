#ifndef AMREX_EBLEVEL_H
#define AMREX_EBLEVEL_H

#include <AMReX_EBLevelGrid.H>
#include <AMReX_EBCellFlag.H>
#include <AMReX_EBFaceFlag.H>

namespace amrex {

class EBLevel
    : public EBLevelGrid
{
public:
    EBLevel ();
    ~EBLevel ();

    EBLevel (const EBLevel& rhs) = default;
    EBLevel (EBLevel&& rhs) noexcept = delete;

    EBLevel& operator= (const EBLevel& rhs) = delete;
    EBLevel& operator= (EBLevel&& rhs) = delete;

    EBLevel (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng);

    void define (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng);

    const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab () const { return *m_cellflags; }
//    const FabArray<EBFaceFlagFab>& getMultiEBFaceFlagFab () const { return *m_faceflags; }
    const EBISLayout&              getEBISL     () const { return *m_ebisl; }


private:
    std::shared_ptr<FabArray<EBCellFlagFab> > m_cellflags;
    std::shared_ptr<FabArray<EBFaceFlagFab> > m_faceflags;
    std::shared_ptr<EBISLayout>               m_ebisl;

    static void setFaceConnection (EBFaceFlag& flag, const FaceIndex& face, const EBISBox& ebis);
    void defineDoit (const BoxArray& ba, const DistributionMapping& dm, const Box& domain, const int ng);
    void setIrregularFaceFab (EBFaceFlagFab & fabs, const EBISBox& ebis, const Box& ccbx);
  
    void setIrregularVolFab (EBCellFlagFab & fab, const EBISBox& ebis, const Box& bx);
};

class MultiFab;
Box getLevelDomain (const MultiFab& mf);

const EBLevel& getEBLevel (const MultiFab& mf);

const FabArray<EBCellFlagFab>& getMultiEBCellFlagFab (const MultiFab& mf);
//const FabArray<EBFaceFlagFab>& getMultiEBFaceFlagFab (const MultiFab& mf);

}

#endif
