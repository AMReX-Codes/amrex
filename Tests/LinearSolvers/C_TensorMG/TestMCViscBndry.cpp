#include <LO_BCTYPES.H>
#include <TestMCViscBndry.H>

void
MCViscBndry::setBndryConds (const BCRec& bc,
			    int          ratio,
			    int          comp)
{
#if BL_SPACEDIM == 2
    BLassert(comp<2*2); // u and v, plus derivs of same
#elif BL_SPACEDIM == 3
    BLassert(comp<3*(3+1)); // u and v, plus derivs of same
#endif

    const REAL* dx = geom.CellSize();
    const BOX& domain = geom.Domain();

    for (OrientationIter fi; fi; ++fi)
    {
	Array<REAL> &bloc = bcloc[fi()];
	Array< Array<BoundCond> >& bctag = bcond[fi()];
	
	int dir = fi().coordDir();
	REAL delta = dx[dir]*ratio;
	int p_bc = (fi().isLow() ? bc.lo(dir): bc.hi(dir));
	
	for (int i = 0; i < boxes().length(); i++)
	{
	    if (domain[fi()] == boxes()[i][fi()] && !geom.isPeriodic(dir))
	    {
		// All physical bc values are located on face
		if (p_bc == EXT_DIR ) {
		    bctag[i][comp] = LO_DIRICHLET;
		    bloc[i] = 0.0;
		} else if (p_bc == FOEXTRAP      ||
			   p_bc == HOEXTRAP      || 
			   p_bc == REFLECT_EVEN)
		{
		    bctag[i][comp] = LO_NEUMANN;
		    bloc[i] = 0.0;
		} else if( p_bc == REFLECT_ODD )
		{
		    bctag[i][comp] = LO_REFLECT_ODD;
		    bloc[i] = 0.0;
		}
	    }
	    else
	    {
		// internal bndry, distance is half of crse
		bctag[i][comp] = LO_DIRICHLET;
		bloc[i] = 0.5*delta;
	    }
	}
    }
}


// *************************************************************************

void
MCViscBndry::setHomogValues()
{
    for (int grd = 0; grd < boxes().length(); grd++)
        for (OrientationIter fi; fi; ++fi)
	    bndry[fi()][grd].setVal(0.);
}
