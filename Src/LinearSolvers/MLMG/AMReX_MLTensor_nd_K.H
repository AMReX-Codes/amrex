#ifndef AMREX_MLTENSOR_ND_K_H_
#define AMREX_MLTENSOR_ND_K_H_

namespace amrex {

  AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
  void mltensor_interpbndry_o1(Box const& box, // area that needs filling
			       int n, // component
			       Array4<Real> const& bcval,
			       Array4<Real const> const& crsebcval,
			       Array4<int const> const&mask, Dim3 const& rr)
{
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    int ic,jc,kc;

    for     (int k = lo.z; k <= hi.z; ++k) {
      kc = ( (rr.z == 0) ? 0 : floor(double(k)/rr.z) );
      for   (int j = lo.y; j <= hi.y; ++j) {
	jc =  ( (rr.y == 0) ? 0 : floor(double(j)/rr.y) );
	for (int i = lo.x; i <= hi.x; i++) {
	    // if mask == 0, cell covered, do nothing
	    if (mask(i,j,k) == BndryData::not_covered) {
	        ic = floor(double(i)/rr.x);
		bcval(i,j,k,n) = crsebcval(ic,jc,kc,n);
	    }
	    // else mask == BndryData::ourside_domain
	    // on physical domain; handled in mltensor_fill
	}
      }
    }
}

}

#endif
