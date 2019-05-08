
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
MLMGBndry::setLOBndryConds (const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& lo,
                            const Vector<Array<LinOpBCType,AMREX_SPACEDIM> >& hi,
                            int ratio, const RealVect& a_loc)
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

        for (int icomp = 0; icomp < nComp(); ++icomp) {
            BCTuple bct;
            setBoxBC(bloc, bct, grd, domain, lo[icomp], hi[icomp], dx, ratio, a_loc);
            for (int idim = 0; idim < 2*AMREX_SPACEDIM; ++idim) {
                bctag[idim][icomp] = bct[idim];
            }
        }
    }
}

void
MLMGBndry::setBoxBC (RealTuple& bloc, BCTuple& bctag, const Box& bx, const Box& domain,
                     const Array<LinOpBCType,AMREX_SPACEDIM>& lo,
                     const Array<LinOpBCType,AMREX_SPACEDIM>& hi,
                     const Real* dx, int ratio, const RealVect& a_loc)
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
                bctag[face] = AMREX_LO_DIRICHLET;
            } else if (linop_bc == LinOpBCType::Neumann) {
                bctag[face] = AMREX_LO_NEUMANN;
            } else if (linop_bc == LinOpBCType::reflect_odd) {
                bctag[face] = AMREX_LO_REFLECT_ODD;
            } else {
                amrex::Abort("MLMGBndry::setBoxBC: Unknown LinOpBCType");
            }
        }
        else
        {
            // Internal bndry.
            bctag[face] = AMREX_LO_DIRICHLET;
            bloc[face]  = ratio > 0 ? 0.5*ratio*dx[dir] : a_loc[dir];
            // If this is next to another same level box, bloc is
            // wrong.  But it doesn't matter, because we also have
            // mask.  It is used only if mask says it is next to
            // coarse cells.
        }
    }
}

}
