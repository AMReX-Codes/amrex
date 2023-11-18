#ifndef AMREX_ML_LIN_OP_SET_H_
#define AMREX_ML_LIN_OP_SET_H_

#include <AMReX_FabArraySet.H>
#include <AMReX_MLLinOp.H>

namespace amrex {

class MLLinOpMFSet
    : public MLLinOpT<MultiFabSet>
{
public:
    //! debug variable
    using MFSet = amrex::MultiFabSet;

    /* Constructors */
    MLLinOpMFSet () = default;
    virtual ~MLLinOpMFSet () = default;

    MLLinOpMFSet (const MLLinOpMFSet&) = delete;
    MLLinOpMFSet (MLLinOpMFSet&&) = delete;
    MLLinOpMFSet& operator= (const MLLinOpMFSet&) = delete;
    MLLinOpMFSet& operator= (MLLinOpMFSet&&) = delete;

    /**** Start public overrides from MLLinOpT ****/
    // void setCoarseFineBC (const MFSet* crse, int crse_ratio) noexcept;
    [[nodiscard]] int getNComp () const override;
    // [[nodiscard]] virtual int getNGrow (int amrlev = 0, int mglev = 0) const { return 0; }
    // void restriction (int amrlev, int cmglev, MFSet& crse, MFSet& fine) const;
    // void interpolation (int amrlev, int fmglev, MFSet& fine, const MFSet& crse) const;
    // void interpAssign (int amrlev, int fmglev, MFSet& fine, MFSet& crse) const;
    // void interpolationAmr (int famrlev, MFSet& fine, const MFSet& crse, const IntVect& nghost) const;
    // void averageDownSolutionRHS (int camrlev, MFSet& crse_sol, MFSet& crse_rhs,
    //                              const MFSet& fine_sol, const MFSet& fine_rhs);
    // virtual void reflux (int crse_amrlev,
    //              MFSet& res, const MFSet& crse_sol, const MFSet& crse_rhs,
    //              MFSet& fine_res, MFSet& fine_sol, const MFSet& fine_rhs) const;
    // RT xdoty (int amrlev, int mglev, const MFSet& x, const MFSet& y, bool local) const;
    // [[nodiscard]] virtual RT normInf (int amrlev, MFSet const& MFSet, bool local) const;
    // virtual void averageDownAndSync (Vector<MFSet>& sol) const;
    // virtual void avgDownResAmr (int clev, MFSet& cres, MFSet const& fres) const;
    // virtual void avgDownResMG (int clev, MFSet& cres, MFSet const& fres) const;
    /**** End public overrides from MLLinOpT ****/

    /**** Start methods inspired by MLCellLinOp and MLNodeLinOp ****/
    // virtual void applyBC (int amrlev, int mglev, MultiFab& phi, BCMode bc_mode, StateMode state_mode,
    //                       bool skip_fillboundary=false) const;
    /**** End methods inspired by MLCellLinOp and MLNodeLinOp ****/

    /**** Start new methods ****/
    int getNSet () const { return m_ixtype_set.size(); }
    /**** End new methods ****/

protected:

    Vector<IndexType> m_ixtype_set;
    int m_ncomp = -1;
    Vector<int> m_ncomp_set;
    Vector<IntVect> m_ngrow_set;

    [[nodiscard]] bool isCellCentered () const noexcept;
    [[nodiscard]] virtual MFSet make (int amrlev, int mglev, IntVect const& ng) const override;
    // [[nodiscard]] virtual MFSet makeAlias (MFSet const& MFSet) const;
    [[nodiscard]] virtual MFSet makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const;
    [[nodiscard]] virtual MFSet makeCoarseAmr (int famrlev, IntVect const& ng) const;
};

using MFSet = amrex::MultiFabSet;

int
MLLinOpMFSet::getNComp () const {
    amrex::Abort("MLLinOpMFSet called on getNComp(). Should be using getNComp(iSet) or getNSet() instead.");
    return m_ncomp;
}

// void
// MLLinOpMFSet::avgDownResMG (int clev, MFSet& cres, MFSet const& fres) const
// {

// #ifdef AMREX_USE_EB
//     if (!fres.isAllRegular()) {
//         amrex::EB_average_down(fres, cres, mg_coarsen_ratio_vec[clev-1]);
//     } else
// #endif
//     {
//         amrex::average_down(fres, cres, mg_coarsen_ratio_vec[clev-1]);
//     }
// }

MFSet
MLLinOpMFSet::make (int amrlev, int mglev, IntVect const& ng) const
{
    auto bxa = m_grids[amrlev][mglev];
    const DistributionMapping& dm = m_dmap[amrlev][mglev];
    auto fact = *m_factory[amrlev][mglev];
    return MFSet(bxa, m_ixtype_set, dm, getNComp(), ng, MFInfo(), fact);
}

// template <typename MFSet>
MFSet
MLLinOpMFSet::makeCoarseMG (int amrlev, int mglev, IntVect const& ng) const
{
    BoxArray cba = m_grids[amrlev][mglev];
    IntVect ratio = (amrlev > 0) ? IntVect(2) : mg_coarsen_ratio_vec[mglev];
    cba.coarsen(ratio);
    return MFSet(cba, m_ixtype_set, m_dmap[amrlev][mglev], getNComp(), ng);
}

// template <typename MFSet>
MFSet
MLLinOpMFSet::makeCoarseAmr (int famrlev, IntVect const& ng) const
{
    BoxArray cba = m_grids[famrlev][0];
    IntVect ratio(AMRRefRatio(famrlev-1));
    cba.coarsen(ratio);
    return MFSet(cba, m_ixtype_set, m_dmap[famrlev][0], getNComp(), ng);
}

}

#endif /* AMREX_ML_LIN_OP_SET_H_ */