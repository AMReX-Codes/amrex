#include "EMParticleContainer.H"
#include "Constants.H"

#include "em_pic_K.H"

using namespace amrex;

void EMParticleContainer::
PushAndDeposeParticles(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                       MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("EMParticleContainer::PushAndDeposeParticles");

    const int lev = 0;

    const auto dxi = Geom(lev).InvCellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();

        ParticleType * pstruct = &(pti.GetArrayOfStructs()[0]);

        auto& attribs = pti.GetAttribs();
        Real *  wp  = attribs[PIdx::w].data();
        Real * uxp  = attribs[PIdx::ux].data();
        Real * uyp  = attribs[PIdx::uy].data();
        Real * uzp  = attribs[PIdx::uz].data();

/*
        Real       * AMREX_RESTRICT Exp  = attribs[PIdx::Ex].data();
        Real       * AMREX_RESTRICT Eyp  = attribs[PIdx::Ey].data();
        Real       * AMREX_RESTRICT Ezp  = attribs[PIdx::Ez].data();
        Real       * AMREX_RESTRICT Bxp  = attribs[PIdx::Bx].data();
        Real       * AMREX_RESTRICT Byp  = attribs[PIdx::By].data();
        Real       * AMREX_RESTRICT Bzp  = attribs[PIdx::Bz].data();
        Real       * AMREX_RESTRICT ginv = attribs[PIdx::ginv].data();
*/

        auto const Exarr = Ex.array(pti);
        auto const Eyarr = Ey.array(pti);
        auto const Ezarr = Ez.array(pti);
        auto const Bxarr = Bx.array(pti);
        auto const Byarr = By.array(pti);
        auto const Bzarr = Bz.array(pti);
        auto       jxarr = jx.array(pti);
        auto       jyarr = jy.array(pti);
        auto       jzarr = jz.array(pti);

        Real q = m_charge;
        Real m = m_mass;
        AMREX_FOR_1D ( np, i,
        {
            amrex::Real Exp;
            amrex::Real Eyp;
            amrex::Real Ezp;
            amrex::Real Bxp;
            amrex::Real Byp;
            amrex::Real Bzp;
            amrex::Real ginv;

            gather_fields(pstruct[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                          Exarr, Eyarr, Ezarr, Bxarr, Byarr, Bzarr, plo, dxi);

            push_momentum_boris(uxp[i], uyp[i], uzp[i], ginv, Exp, Eyp, Ezp,
                                Bxp, Byp, Bzp, q, m, dt);

            push_position_boris(pstruct[i], uxp[i], uyp[i], uzp[i], ginv, dt);

            deposit_current(jxarr, jyarr, jzarr, pstruct[i], uxp[i], uyp[i], uzp[i],
                            ginv, wp[i], q, dt, plo, dxi);
        });       
    }
}

void EMParticleContainer::
PushParticleMomenta(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                    const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz, Real dt)
{   
    BL_PROFILE("EMParticleContainer::PushParticleMomenta");

    const int lev = 0;
    
    const auto dxi = Geom(lev).InvCellSizeArray();
    const auto plo = Geom(lev).ProbLoArray();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();
        
        ParticleType const* AMREX_RESTRICT pstruct = &(pti.GetArrayOfStructs()[0]);

        auto& attribs = pti.GetAttribs();
        Real * AMREX_RESTRICT uxp  = attribs[PIdx::ux].data();
        Real * AMREX_RESTRICT uyp  = attribs[PIdx::uy].data();
        Real * AMREX_RESTRICT uzp  = attribs[PIdx::uz].data();
        /*
        Real * AMREX_RESTRICT Exp  = attribs[PIdx::Ex].data();
        Real * AMREX_RESTRICT Eyp  = attribs[PIdx::Ey].data();
        Real * AMREX_RESTRICT Ezp  = attribs[PIdx::Ez].data();
        Real * AMREX_RESTRICT Bxp  = attribs[PIdx::Bx].data();
        Real * AMREX_RESTRICT Byp  = attribs[PIdx::By].data();
        Real * AMREX_RESTRICT Bzp  = attribs[PIdx::Bz].data();
        Real * AMREX_RESTRICT ginv = attribs[PIdx::ginv].data();
        */

        auto const Exarr = Ex.array(pti);
        auto const Eyarr = Ey.array(pti);
        auto const Ezarr = Ez.array(pti);
        auto const Bxarr = Bx.array(pti);
        auto const Byarr = By.array(pti);
        auto const Bzarr = Bz.array(pti);

        Real q = m_charge;
        Real m = m_mass;

        AMREX_PARALLEL_FOR_1D ( np, i,
        {
            amrex::Real Exp;
            amrex::Real Eyp;
            amrex::Real Ezp;
            amrex::Real Bxp;
            amrex::Real Byp;
            amrex::Real Bzp;
            amrex::Real ginv;
            gather_fields(pstruct[i], Exp, Eyp, Ezp, Bxp, Byp, Bzp,
                          Exarr, Eyarr, Ezarr, Bxarr, Byarr, Bzarr, plo, dxi);

            push_momentum_boris(uxp[i], uyp[i], uzp[i], ginv, Exp, Eyp, Ezp,
                                Bxp, Byp, Bzp, q, m, dt);
        });
    }
}

