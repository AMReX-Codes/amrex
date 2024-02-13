#ifndef AMREX_GMRES_MLMG_H_
#define AMREX_GMRES_MLMG_H_
#include <AMReX_Config.H>

#include <AMReX_MLMG.H>
#include <utility>

namespace amrex {

//! Wrapping MLMG as a matrix operator for GMRES
template <typename M>
class GMRESMLMGT
{
public:
    using MF = typename M::MFType; // typically MultiFab
    using RT = typename MF::value_type; // double or float

    explicit GMRESMLMGT (M& mlmg);

    //! Make MultiFab without ghost cells
    MF makeVecRHS () const;

    //! Make MultiFab with ghost cells and set ghost cells to zero
    MF makeVecLHS () const;

    RT norm2 (MF const& mf) const;

    static void scale (MF& mf, RT scale_factor);

    RT dotProduct (MF const& mf1, MF const& mf2) const;

    //! lhs = value
    static void setVal (MF& lhs, RT value);

    //! lhs = rhs
    static void assign (MF& lhs, MF const& rhs);

    //! lhs += a*rhs
    static void increment (MF& lhs, MF const& rhs, RT a);

    //! lhs = a*rhs_a + b*rhs_b
    static void linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b);

    //! lhs = L(rhs)
    void apply (MF& lhs, MF const& rhs) const;

    void precond (MF& lhs, MF const& rhs) const;

    bool usePrecond (bool new_flag) { return std::exchange(m_use_precond, new_flag); }

private:
    M& m_mlmg;
    MLLinOpT<MF>& m_linop;
    bool m_use_precond = false;
};

template <typename M>
GMRESMLMGT<M>::GMRESMLMGT (M& mlmg)
    : m_mlmg(mlmg), m_linop(mlmg.getLinOp())
{
    m_mlmg.prepareLinOp();
}

template <typename M>
auto GMRESMLMGT<M>::makeVecRHS () const -> MF
{
    return m_linop.make(0, 0, IntVect(0));
}

template <typename M>
auto GMRESMLMGT<M>::makeVecLHS () const -> MF
{
    auto mf = m_linop.make(0, 0, IntVect(1));
    mf.setBndry(0);
    return mf;
}

template <typename M>
auto GMRESMLMGT<M>::norm2 (MF const& mf) const -> RT
{
    auto r = m_linop.xdoty(0, 0, mf, mf, false);
    return std::sqrt(r);
}

template <typename M>
void GMRESMLMGT<M>::scale (MF& mf, RT scale_factor)
{
    mf.mult(scale_factor, 0, mf.nComp());
}

template <typename M>
auto GMRESMLMGT<M>::dotProduct (MF const& mf1, MF const& mf2) const -> RT
{
    return m_linop.xdoty(0, 0, mf1, mf2, false);
}

template <typename M>
void GMRESMLMGT<M>::setVal (MF& lhs, RT value)
{
    lhs.setVal(value);
}

template <typename M>
void GMRESMLMGT<M>::assign (MF& lhs, MF const& rhs)
{
    MF::Copy(lhs, rhs, 0, 0, lhs.nComp(), IntVect(0));
}

template <typename M>
void GMRESMLMGT<M>::increment (MF& lhs, MF const& rhs, RT a)
{
    MF::Saxpy(lhs, a, rhs, 0, 0, lhs.nComp(), IntVect(0));
}

template <typename M>
void GMRESMLMGT<M>::linComb (MF& lhs, RT a, MF const& rhs_a, RT b, MF const& rhs_b)
{
    MF::LinComb(lhs, a, rhs_a, 0, b, rhs_b, 0, 0, lhs.nComp(), IntVect(0));
}

template <typename M>
void GMRESMLMGT<M>::apply (MF& lhs, MF const& rhs) const
{
    m_linop.apply(0, 0, lhs, const_cast<MF&>(rhs),
                  MLLinOpT<MF>::BCMode::Homogeneous,
                  MLLinOpT<MF>::StateMode::Correction);
}

template <typename M>
void GMRESMLMGT<M>::precond (MF& lhs, MF const& rhs) const
{
    if (m_use_precond) {
        // for now, let's just do some smoothing
        lhs.setVal(RT(0.0));
        for (int m = 0; m < 4; ++m) {
            m_linop.smooth(0, 0, lhs, rhs, (m==0) ? true : false);
        }
    } else {
        amrex::Copy(lhs, rhs, 0, 0, lhs.nComp(), IntVect(0));
    }
}

using GMRESMLMG = GMRESMLMGT<MLMG>;

}

#endif
