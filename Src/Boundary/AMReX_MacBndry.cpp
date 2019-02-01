
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MacBndry.H>

namespace amrex {

MacBndry::MacBndry ()
    :
    InterpBndryData()
{
    amrex::Abort("*** Calling default constructor for MacBndry()");
}

MacBndry::MacBndry (const BoxArray& _grids,
		    const DistributionMapping& _dmap,
                    int             _ncomp,
                    const Geometry& _geom)
    :
    InterpBndryData(_grids,_dmap,_ncomp,_geom)
{}

MacBndry::~MacBndry () {}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         int          ratio)
{
    const IntVect& ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
MacBndry::setBndryConds (const BCRec&   phys_bc,
                         const IntVect& ratio,
			 int            comp)
{
    m_phys_bc = phys_bc;

    //
    // ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    // DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const BoxArray& ba     = boxes();
    const Real*     dx     = geom.CellSize();
    const Box&      domain = geom.Domain();
    //
    // We note that all orientations of the FabSets have the same distribution.
    // We'll use the low 0 side as the model.
    //
    for (FabSetIter fsi(bndry[Orientation(0,Orientation::low)]); fsi.isValid(); ++fsi)
    {
        const int                  i     = fsi.index();
        const Box&                 grd   = ba[i];
        RealTuple&                 bloc  = bcloc[fsi];
        Vector< Vector<BoundCond> >& bctag = bcond[fsi];

        for (OrientationIter fi; fi; ++fi)
        {
            const Orientation face  = fi();
            const int         dir   = face.coordDir();

            if (domain[face] == grd[face] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                const int p_bc  = (face.isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

                bctag[face][comp] = (p_bc == PhysBCType::outflow) 
                    ? AMREX_LO_DIRICHLET : AMREX_LO_NEUMANN;
                bloc[face]        = 0;
            }
            else
            {
                //
                // Internal bndry.
                //
                const Real delta = dx[dir]*ratio[dir];

                bctag[face][comp] = AMREX_LO_DIRICHLET;
		bloc[face]        = 0.5*delta;
            }
        }
    }
}

void
MacBndry::setHomogValues (const BCRec&   bc,
                          const IntVect& ratio)
{
    setBndryConds(bc, ratio);
    setVal(0.0);
}

}
