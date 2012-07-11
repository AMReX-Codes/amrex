
#include <winstd.H>
#include <LO_BCTYPES.H>
#include <MacBndry.H>

MacBndry::MacBndry ()
    :
    InterpBndryData()
{}

MacBndry::MacBndry (const BoxArray& _grids,
                    int             _ncomp,
                    const Geometry& _geom)
    :
    InterpBndryData(_grids,_ncomp,_geom)
{}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         int          ratio)
{
    IntVect ratio_vect = ratio * IntVect::TheUnitVector();
    setBndryConds(phys_bc, ratio_vect);
}

void
MacBndry::setBndryConds (const BCRec& phys_bc,
                         IntVect&     ratio,
			 int          comp)
{
    //
    // ALL BCLOC VALUES ARE NOW DEFINED AS A LENGTH IN PHYSICAL
    // DIMENSIONS *RELATIVE* TO THE FACE, NOT IN ABSOLUTE PHYSICAL SPACE
    //
    const BoxArray& grids  = boxes();
    const Real*     dx     = geom.CellSize();
    const Box&      domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
        Array<Real>& bloc = bcloc[fi()];

        const int  dir   = fi().coordDir();
        const Real delta = dx[dir]*ratio[dir];
        const int  p_bc  = (fi().isLow() ? phys_bc.lo(dir) : phys_bc.hi(dir));

        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            const int i = fsi.index();

            const Box& grd = grids[i];

            Array<BoundCond>& bctag = bcond[fi()][i];

            if (domain[fi()] == grd[fi()] && !geom.isPeriodic(dir))
            {
                //
                // All physical bc values are located on face.
                //
                bctag[comp] = (p_bc == Outflow) ? LO_DIRICHLET : LO_NEUMANN;
                bloc[i]        = 0;
            }
            else
            {
                //
                // Internal bndry.
                //
                bctag[comp] = LO_DIRICHLET;
		bloc[i] = 0.5*delta;
            }
        }
    }
}

void
MacBndry::setHomogValues (const BCRec& bc,
                          /*const*/ IntVect& ratio)
{
    setBndryConds(bc, ratio);
 
    for (OrientationIter fi; fi; ++fi)
    {
        for (FabSetIter fsi(bndry[fi()]); fsi.isValid(); ++fsi)
        {
            bndry[fi()][fsi].setVal(0);
        }
    }
}

