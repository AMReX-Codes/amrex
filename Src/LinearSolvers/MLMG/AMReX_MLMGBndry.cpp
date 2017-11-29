
#include <AMReX_MLMGBndry.H>

namespace amrex {

MLMGBndry::MLMGBndry (const BoxArray& _grids,
                      const DistributionMapping& _dmap,
                      int             _ncomp,
                      const Geometry& _geom)
    : InterpBndryData(_grids,_dmap,_ncomp,_geom)
{}

MLMGBndry::~MLMGBndry () {}

void
MLMGBndry::setLOBndryConds (const std::array<LinOpBCType,AMREX_SPACEDIM>& lo,
                            const std::array<LinOpBCType,AMREX_SPACEDIM>& hi,
                            int ratio)
{
    const BoxArray& ba     = boxes();
    const Real*     dx     = geom.CellSize();
    const Box&      domain = geom.Domain();

#ifdef _OPENMP
#pragma omp parallel
#endif
    for (FabSetIter fsi(bndry[Orientation(0,Orientation::low)]); fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        const Box&                 grd   = ba[i];
        RealTuple&                 bloc  = bcloc[fsi];
        Vector< Vector<BoundCond> >& bctag = bcond[fsi];

        BCTuple bct;
        setBoxBC(bloc, bct, grd, domain, lo, hi, dx, ratio);

        const int comp = 0;
        for (int idim = 0; idim < 2*AMREX_SPACEDIM; ++idim) {
            bctag[idim][comp] = bct[idim];
        }
    }
}

void
MLMGBndry::setBoxBC (RealTuple& bloc, BCTuple& bctag, const Box& bx, const Box& domain,
                     const std::array<LinOpBCType,AMREX_SPACEDIM>& lo,
                     const std::array<LinOpBCType,AMREX_SPACEDIM>& hi,
                     const Real* dx, int ratio)
{
    for (OrientationIter fi; fi; ++fi)
    {
        const Orientation face = fi();
        const int         dir  = face.coordDir();
        
        if (domain[face] == bx[face] && !Geometry::isPeriodic(dir))
        {
            // All physical bc values are located on face.
            bloc[face]  = 0;
            const auto linop_bc  = face.isLow() ? lo[dir] : hi[dir];
            if (linop_bc == LinOpBCType::Dirichlet) {
                bctag[face] = LO_DIRICHLET;
            } else if (linop_bc == LinOpBCType::Neumann) {
                bctag[face] = LO_NEUMANN;
            } else if (linop_bc == LinOpBCType::reflect_odd) {
                bctag[face] = LO_REFLECT_ODD;
            } else {
                amrex::Abort("MLMGBndry::setBoxBC: Unknown LinOpBCType");
            }
        }
        else
        {
            // Internal bndry.
            bctag[face] = LO_DIRICHLET;
            bloc[face]  = 0.5*ratio*dx[dir];
            // If this is next to another same level box, bloc is
            // wrong.  But it doesn't matter, because we also have
            // mask.  It is used only if mask says it is next to
            // coarse cells.
        }
    }
}

}
