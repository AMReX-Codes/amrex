#include "EMParticleContainer.H"
#include "Constants.H"

#include "em_pic_F.H"

using namespace amrex;

void EMParticleContainer::
PushAndDeposeParticles(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                       const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz,
                       MultiFab& jx, MultiFab& jy, MultiFab& jz, Real dt)
{
    BL_PROFILE("EMParticleContainer::PushAndDeposeParticles");

    const int lev = 0;

    const Real* dx  = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();

        auto& structs = pti.GetArrayOfStructs();
        auto& attribs = pti.GetAttribs();

        auto&  wp  = attribs[PIdx::w];
        auto& uxp  = attribs[PIdx::ux];
        auto& uyp  = attribs[PIdx::uy];
        auto& uzp  = attribs[PIdx::uz];
        auto& Exp  = attribs[PIdx::Ex];
        auto& Eyp  = attribs[PIdx::Ey];
        auto& Ezp  = attribs[PIdx::Ez];
        auto& Bxp  = attribs[PIdx::Bx];
        auto& Byp  = attribs[PIdx::By];
        auto& Bzp  = attribs[PIdx::Bz];
        auto& ginv = attribs[PIdx::ginv];

	gather_magnetic_field(np, structs.dataPtr(),
                                    Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                    BL_TO_FORTRAN_3D(Bx[pti]),
                                    BL_TO_FORTRAN_3D(By[pti]),
                                    BL_TO_FORTRAN_3D(Bz[pti]),
                                    plo, dx);

	gather_electric_field(np, structs.dataPtr(),
                                    Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                    BL_TO_FORTRAN_3D(Ex[pti]),
                                    BL_TO_FORTRAN_3D(Ey[pti]),
                                    BL_TO_FORTRAN_3D(Ez[pti]),
                                    plo, dx);

        push_momentum_boris(np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), ginv.dataPtr(),
                                  Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                  Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                  m_charge, m_mass, dt);

        push_position_boris(np, structs.dataPtr(),
                                  uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), ginv.dataPtr(), dt);

        deposit_current(BL_TO_FORTRAN_3D(jx[pti]),
                              BL_TO_FORTRAN_3D(jy[pti]),
                              BL_TO_FORTRAN_3D(jz[pti]),
                              np, structs.dataPtr(),
                              uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(),
                              ginv.dataPtr(), wp.dataPtr(),
                              m_charge, plo, dt, dx);
    }
}

void EMParticleContainer::
PushParticleMomenta(const MultiFab& Ex, const MultiFab& Ey, const MultiFab& Ez,
                    const MultiFab& Bx, const MultiFab& By, const MultiFab& Bz, Real dt)
{
    BL_PROFILE("EMParticleContainer::PushParticleMomenta");

    const int lev = 0;

    const Real* dx  = Geom(lev).CellSize();
    const Real* plo = Geom(lev).ProbLo();

    for (EMParIter pti(*this, lev); pti.isValid(); ++pti)
    {
        const int np  = pti.numParticles();

        auto& structs = pti.GetArrayOfStructs();
        auto& attribs = pti.GetAttribs();

        auto& uxp  = attribs[PIdx::ux];
        auto& uyp  = attribs[PIdx::uy];
        auto& uzp  = attribs[PIdx::uz];
        auto& Exp  = attribs[PIdx::Ex];
        auto& Eyp  = attribs[PIdx::Ey];
        auto& Ezp  = attribs[PIdx::Ez];
        auto& Bxp  = attribs[PIdx::Bx];
        auto& Byp  = attribs[PIdx::By];
        auto& Bzp  = attribs[PIdx::Bz];
        auto& ginv = attribs[PIdx::ginv];

	gather_magnetic_field(np, structs.dataPtr(),
                                    Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                    BL_TO_FORTRAN_3D(Bx[pti]),
                                    BL_TO_FORTRAN_3D(By[pti]),
                                    BL_TO_FORTRAN_3D(Bz[pti]),
                                    plo, dx);

	gather_electric_field(np, structs.dataPtr(),
                                    Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                    BL_TO_FORTRAN_3D(Ex[pti]),
                                    BL_TO_FORTRAN_3D(Ey[pti]),
                                    BL_TO_FORTRAN_3D(Ez[pti]),
                                    plo, dx);

        push_momentum_boris(np, uxp.dataPtr(), uyp.dataPtr(), uzp.dataPtr(), ginv.dataPtr(),
                                  Exp.dataPtr(), Eyp.dataPtr(), Ezp.dataPtr(),
                                  Bxp.dataPtr(), Byp.dataPtr(), Bzp.dataPtr(),
                                  m_charge, m_mass, dt);
    }
}
