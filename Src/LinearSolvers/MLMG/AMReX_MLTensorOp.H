#ifndef AMREX_ML_TENSOR_OP_H_
#define AMREX_ML_TENSOR_OP_H_

#include <AMReX_MLCellABecLap.H>
#include <AMReX_Array.H>

namespace amrex {

// Tensor solver for high Reynolds flows with near constant shear viscosity.
// 
// a v + [-div dot (eta grad v) + (grad eta) (div dot v) - (grad eta) dot (grad v)^T] = rhs
//
// Note that what's inside [] is div dot tau + div dot ((kappa+eta/3) I (div dot v)),
// where tau is viscous tensor, and eta and kappa are shear and buld viscosity, respectively.

class MLTensorOp
    : public MLCellABecLap
{
public:

    MLTensorOp () {}
    MLTensorOp (const Vector<Geometry>& a_geom,
                const Vector<BoxArray>& a_grids,
                const Vector<DistributionMapping>& a_dmap,
                const LPInfo& a_info = LPInfo(),
                const Vector<FabFactory<FArrayBox> const*>& a_factory = {});
    virtual ~MLTensorOp ();

    MLTensorOp (const MLTensorOp&) = delete;
    MLTensorOp (MLTensorOp&&) = delete;
    MLTensorOp& operator= (const MLTensorOp&) = delete;
    MLTensorOp& operator= (MLTensorOp&&) = delete;

    void define (const Vector<Geometry>& a_geom,
                 const Vector<BoxArray>& a_grids,
                 const Vector<DistributionMapping>& a_dmap,
                 const LPInfo& a_info = LPInfo(),
                 const Vector<FabFactory<FArrayBox> const*>& a_factory = {});

    void setACoeffs (int amrlev, const MultiFab& alpha);
    void setBCoeffs (int amrlev, const Array<MultiFab const*,AMREX_SPACEDIM>& beta);

    virtual void prepareForSolve () final override;
    virtual bool isSingular (int armlev) const final override { return false; }
    virtual bool isBottomSingular () const final override { return false; }
    virtual void Fapply (int amrlev, int mglev, MultiFab& out, const MultiFab& in) const final override;
    virtual void Fsmooth (int amrlev, int mglev, MultiFab& sol, const MultiFab& rhs, int redblack) const final override;
    virtual void FFlux (int amrlev, const MFIter& mfi,
                        const Array<FArrayBox*,AMREX_SPACEDIM>& flux,
                        const FArrayBox& sol, Location /* loc */,
                        const int face_only=0) const final override;

    virtual void normalize (int amrlev, int mglev, MultiFab& mf) const final override;

    virtual Real getAScalar () const final override { return 1.0; }
    virtual Real getBScalar () const final override { return 1.0; }
    virtual MultiFab const* getACoeffs (int amrlev, int mglev) const final override
        { return &(m_a_coeffs[amrlev][mglev]); }
    virtual Array<MultiFab const*,AMREX_SPACEDIM> getBCoeffs (int amrlev, int mglev) const final override
        { return amrex::GetArrOfConstPtrs(m_b_coeffs[amrlev][mglev]); }

    virtual std::unique_ptr<MLLinOp> makeNLinOp (int grid_size) const final override {
        amrex::Abort("MLTensorOp::makeNLinOp: Not implmented");
        return std::unique_ptr<MLLinOp>{};
    }

    void averageDownCoeffsSameAmrLevel (Vector<MultiFab>& a,
                                        Vector<Array<MultiFab,AMREX_SPACEDIM> >& b);
    void averageDownCoeffs ();
    void averageDownCoeffsToCoarseAmrLevel (int flev);

private:

    Vector<Vector<MultiFab> > m_a_coeffs;
    Vector<Vector<Array<MultiFab,AMREX_SPACEDIM> > > m_b_coeffs;
};

}

#endif
