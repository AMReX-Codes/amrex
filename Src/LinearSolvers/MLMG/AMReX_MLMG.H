#ifndef AMREX_ML_MG_H_
#define AMREX_ML_MG_H_

#include <AMReX_MLLinOp.H>

namespace amrex {

class MLMG
{
public:

    MLMG (MLLinOp& a_lp);
    ~MLMG ();

    void solve (const Vector<MultiFab*>& sol, const Vector<MultiFab*>& rhs,
                Real a_tol_real, Real a_tol_abs);

private:
    
    MLLinOp& m_lp;

};

}

#endif
