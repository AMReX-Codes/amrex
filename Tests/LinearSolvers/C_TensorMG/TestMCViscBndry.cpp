#include <LO_BCTYPES.H>
#include <TestMCViscBndry.H>

void MCViscBndry::setBndryConds(const Array<BCRec>& bcarray,
			      const Geometry& geom, int ratio)
{
    int ncomp = bcond[0][0].length();
#if BL_SPACEDIM == 2
    assert(ncomp==2*2); // u and v, plus derivs of same
    assert(bcarray.length()==2*2);
#elif BL_SPACEDIM == 3
    assert(ncomp==3*(3+1)); // u and v, plus derivs of same
    assert(bcarray.length()==3*(3+1));
#endif

    const BoxArray& grids = boxes();
    int ngrds = grids.length();
    const REAL* dx = geom.CellSize();
    const BOX& domain = geom.Domain();
    const RealBox& prob_domain = geom.ProbDomain();

    for (OrientationIter fi; fi; ++fi) {
	Orientation face(fi());
	Array<REAL> &bloc = bcloc[face];

	int dir = face.coordDir();
	REAL delta = dx[dir]*ratio;
	for( int icomp=0; icomp<ncomp; icomp++){
	  int p_bc = (face.isLow() ? 
		      bcarray[icomp].lo(dir):bcarray[icomp].hi(dir));

	  for (int i = 0; i < ngrds; i++) {
	    const BOX& grd = grids[i];
	    // bctag is bc type (with array info on orientation,grid,comp)gone
	    BoundCond  &bctag = bcond[face][i][icomp];

	    if (domain[face] == grd[face]) {
	      // All physical bc values are located on face
	      if (p_bc == EXT_DIR ) {
		bctag = LO_DIRICHLET;
		bloc[i] = 0.0; // on face, distance to face = 0
	      } else if (p_bc == EXTRAP || p_bc == HOEXTRAP || 
			 p_bc == REFLECT_EVEN) {
		bctag = LO_NEUMANN;
		bloc[i] = 0.0; // on face, distance to face = 0
	      } else if( p_bc == REFLECT_ODD ){
		bctag = LO_REFLECT_ODD;
		bloc[i] = 0.0; // on face, distance to face = 0
	      }
	    } else {
	      // internal bndry
	      bctag = LO_DIRICHLET;
	      bloc[i] = 0.5*delta; // internal, distance is half of crse
	    }
	  }
	}
    }
}


// *************************************************************************

void
MCViscBndry::setHomogValues()
{
    int ngrd = grids.length();
    for (int grd = 0; grd < ngrd; grd++) {
        const BOX& bx = grids[grd];
        for (OrientationIter fi; fi; ++fi) {
            Orientation face(fi());
            FARRAYBOX& bnd_fab = bndry[face][grd];
            bnd_fab.setVal(0.);
        }
    }
}
