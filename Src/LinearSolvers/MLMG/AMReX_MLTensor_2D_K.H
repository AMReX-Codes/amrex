#ifndef AMREX_MLTENSOR_2D_K_H_
#define AMREX_MLTENSOR_2D_K_H_

#include <AMReX_MLLinOp_K.H>
#include <AMReX_MLTensor_nd_K.H>

namespace amrex {
  
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_fill_corners (int icorner, Box const& vbox, // vbox: the valid box
                            Array4<Real> const& vel,
                            Array4<int const> const& mxlo,
                            Array4<int const> const& mylo,
                            Array4<int const> const& mxhi,
                            Array4<int const> const& myhi,
                            Array4<Real> const& bcvalxlo,
                            Array4<Real> const& bcvalylo,
                            Array4<Real> const& bcvalxhi,
                            Array4<Real> const& bcvalyhi,
                            Array4<Real const> const& crsebcvalxlo,
                            Array4<Real const> const& crsebcvalylo,
                            Array4<Real const> const& crsebcvalxhi,
                            Array4<Real const> const& crsebcvalyhi,
                            GpuArray<BoundCond,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bct,
                            GpuArray<Real,2*AMREX_SPACEDIM*AMREX_SPACEDIM> const& bcl,
                            int inhomog, int maxorder, int use_crsedata,
                            GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
			    Dim3 const& rr,
                            Box const& domain) noexcept
{
    constexpr int oxlo = 0;
    constexpr int oylo = 1;
    constexpr int oxhi = 2;
    constexpr int oyhi = 3;
    constexpr int xdir = 0;
    constexpr int ydir = 1;
    const auto blen = amrex::length(vbox);
    const auto vlo  = amrex::lbound(vbox);
    const auto vhi  = amrex::ubound(vbox);
    const auto dlo  = amrex::lbound(domain);
    const auto dhi  = amrex::ubound(domain);

    //fixme
    //static int count=0; count++;
    //
    for (int icomp = 0; icomp < AMREX_SPACEDIM; ++icomp) {
        switch (icorner) {
        case 0: {
            // xlo & ylo
            Box bx = amrex::adjCellLo(amrex::adjCellLo(vbox,xdir,1),ydir,1);
	    //const auto blo  = amrex::lbound(bx);
            if (vlo.x == dlo.x) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalxlo, crsebcvalxlo, mxlo, rr);

                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            }
            if (vlo.y == dlo.y) {
	        int offset = AMREX_SPACEDIM * oylo;
		if ( use_crsedata )
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalylo, crsebcvalylo, mylo, rr);
		
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
 	    // if (vbox.contains({0,0})){
	    //   int i = 0, j=0;
	    //   Print()<<"count "<<count<<": vbx "<<vbox<<"\n";
	    //   Print()<<"count "<<count<<": bcval("<<i-1<<","<<j-1<<","<<icomp<<") "<<bcvalylo(i-1,j-1,0,icomp)<<" "<<bcvalxlo(i-1,j-1,0,icomp)<<" "<<vel(i-1,j-1,0,icomp)<<"\n";
	    //   Print()<<"count "<<count<<": bcval("<<icomp<<") "<<bcvalylo(i,j-1,0,icomp)<<" "<<bcvalxlo(i-1,j,0,icomp)<<"\n";
	    //   //amrex::Abort();
	    // }
            }
            break;
        }
        case 1: {
            // xhi & ylo
            Box bx = amrex::adjCellLo(amrex::adjCellHi(vbox,xdir,1),ydir,1);
            if (vhi.x == dhi.x) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalxhi, crsebcvalxhi, mxhi, rr);

                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            }
            if (vlo.y == dlo.y) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalylo, crsebcvalylo, mylo, rr);

                int offset = AMREX_SPACEDIM * oylo;
                mllinop_apply_bc_y(Orientation::low, bx, blen.y,
                                   vel, mylo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalylo, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        case 2: {
            // xlo & yhi
            Box bx = amrex::adjCellHi(amrex::adjCellLo(vbox,xdir,1),ydir,1);
            if (vlo.x == dlo.x) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalxlo, crsebcvalxlo, mxlo, rr);

                int offset = AMREX_SPACEDIM * oxlo;
                mllinop_apply_bc_x(Orientation::low, bx, blen.x,
                                   vel, mxlo, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxlo, maxorder, dxinv[xdir], inhomog, icomp);
            }
            if (vhi.y == dhi.y) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalyhi, crsebcvalyhi, myhi, rr);

                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            }

	    //Print()<<ParallelDescriptor::MyProc()<<": call "<<count<<"\n";

	    // if (vbox.contains({128,392})){
	    // 	  Print()<<"bx "<<bx<<"\n";
	    // 	  Print()<<"dlo "<<dlo<<"\n";
	    // 	  Print()<<"bcvalxlo "<<bcvalxlo(127,391,0,icomp);
	    // // }
	    // if (vbox.contains({16,47})){
	    //   Print()<<"bx "<<bx<<"\n";
	    //   Print()<<"dlo "<<dlo<<"\n";
	    //   Print()<<"bcvalyhi "<<bcvalyhi(15,48,0,icomp)<<" "<<bcvalyhi(16,48,0,icomp)<<"\n";
	    // }
	    // if (vbox.contains({15,47})){
	    //   Print()<<"vbx "<<vbox<<"\n";
	    //   Print()<<"bcvalxlo "<<bcvalxlo(-1,31,0,icomp)<<" "<<bcvalxlo(-1,32,0,icomp)<<" "<<bcvalxlo(-1,47,0,icomp)<<" "<<bcvalxlo(-1,48,0,icomp)<<"\n";
	    //   Print()<<"bcvalylo "<<bcvalylo(-1,31,0,icomp)<<" "<<bcvalylo(0,31,0,icomp)<<" "<<bcvalylo(15,31,0,icomp)<<" "<<bcvalylo(16,31,0,icomp)<<"\n";
	    //   Print()<<"bcvalxhi "<<bcvalxhi(16,31,0,icomp)<<" "<<bcvalxhi(16,32,0,icomp)<<" "<<bcvalxhi(16,47,0,icomp)<<" "<<bcvalxhi(16,48,0,icomp)<<"\n";
	    //   Print()<<"bcvalyhi "<<bcvalyhi(-1,48,0,icomp)<<" "<<bcvalyhi(0,48,0,icomp)<<" "<<bcvalyhi(15,48,0,icomp)<<" "<<bcvalyhi(16,48,0,icomp)<<"\n";
	    // }


            break;
        }
        case 3: {
            // xhi & yhi
            Box bx = amrex::adjCellHi(amrex::adjCellHi(vbox,xdir,1),ydir,1);
            if (vhi.x == dhi.x) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalxhi, crsebcvalxhi, mxhi, rr);
	      
                int offset = AMREX_SPACEDIM * oxhi;
                mllinop_apply_bc_x(Orientation::high, bx, blen.x,
                                   vel, mxhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalxhi, maxorder, dxinv[xdir], inhomog, icomp);
            }
            if (vhi.y == dhi.y) {
	        if (use_crsedata)
		  // fill corners of bndry from coarse data using first order interpolation 
		  mltensor_interpbndry_o1(bx, icomp, bcvalyhi, crsebcvalyhi, myhi, rr);

                int offset = AMREX_SPACEDIM * oyhi;
                mllinop_apply_bc_y(Orientation::high, bx, blen.y,
                                   vel, myhi, bct[offset+icomp], bcl[offset+icomp],
                                   bcvalyhi, maxorder, dxinv[ydir], inhomog, icomp);
            }
            break;
        }
        default: {}
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fx (Box const& box, Array4<Real> const& fx,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etax,
                              Array4<Real const> const& kapx,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dyi = dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudy = (vel(i,j+1,0,0)+vel(i-1,j+1,0,0)-vel(i,j-1,0,0)-vel(i-1,j-1,0,0))*(0.25*dyi);
            Real dvdy = (vel(i,j+1,0,1)+vel(i-1,j+1,0,1)-vel(i,j-1,0,1)-vel(i-1,j-1,0,1))*(0.25*dyi);
            //Real divu = dvdy;
	    //fixme
	    Real divu = 0.;
            Real xif = kapx(i,j,0);
            Real mun = 0.75*(etax(i,j,0,0)-xif);  // restore the original eta
            Real mut =       etax(i,j,0,1);
            fx(i,j,0,0) = -mun*(-twoThirds*divu) - xif*divu;
            fx(i,j,0,1) = -mut*dudy;
        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms_fy (Box const& box, Array4<Real> const& fy,
                              Array4<Real const> const& vel,
                              Array4<Real const> const& etay,
                              Array4<Real const> const& kapy,
                              GpuArray<Real,AMREX_SPACEDIM> const& dxinv) noexcept
{
    const Real dxi = dxinv[0];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);
    constexpr Real twoThirds = 2./3.;

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Real dudx = (vel(i+1,j,0,0)+vel(i+1,j-1,0,0)-vel(i-1,j,0,0)-vel(i-1,j-1,0,0))*(0.25*dxi);
            Real dvdx = (vel(i+1,j,0,1)+vel(i+1,j-1,0,1)-vel(i-1,j,0,1)-vel(i-1,j-1,0,1))*(0.25*dxi);
            //Real divu = dudx;
	    //fixme
	    Real divu = 0.;
            Real xif = kapy(i,j,0);
            Real mun = 0.75*(etay(i,j,0,1)-xif);  // restore the original eta
            Real mut =       etay(i,j,0,0);
            fy(i,j,0,0) = -mut*dvdx;
            fy(i,j,0,1) = -mun*(-twoThirds*divu) - xif*divu;
	    //fixme
	    // IntVect mypt={0,0};
	    // if (box.contains(mypt)) {
	    // 	if (i==0 && j==0){
	    // 	  Print()<<" i,j: "<<i<< ","<<j<<"\n";
	    // 	  Print()<<"vel(-1,-1) "<<vel(i-1,j-1,0,0)<<" "<<vel(i-1,j-1,0,1)<<"\n";
	    // 	  Print()<<"vel(-1,0) "<<vel(i-1,j,0,0)<<" "<<vel(i-1,j,0,1)<<"\n";
	    // 	  Print()<<"vel(0,0) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n";
	    // 	  Print()<<"vel(0,-1) "<<vel(i,j-1,0,0)<<" "<<vel(i,j-1,0,1)<<"\n";
	    // 	}
	    // 	if (i==0 && j==15){
	    // 	  Print()<<" i,j: "<<i<< ","<<j<<"\n";
	    // 	  Print()<<"vel(-1,-1) "<<vel(i-1,j-1,0,0)<<" "<<vel(i-1,j-1,0,1)<<"\n";
	    // 	  Print()<<"vel(-1,0) "<<vel(i-1,j,0,0)<<" "<<vel(i-1,j,0,1)<<"\n";
	    // 	  Print()<<"vel(0,0) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n";
	    // 	  Print()<<"vel(0,-1) "<<vel(i,j-1,0,0)<<" "<<vel(i,j-1,0,1)<<"\n";
	    // 	}
	    //   }
	    //  mypt={0,47};
	    // if (box.contains(mypt)) {
	    	// if (i==0 && j==32){
	    	//   Print()<<" i,j: "<<i<< ","<<j<<"\n";
	    	//   Print()<<"vel(-1,-1) "<<vel(i-1,j-1,0,0)<<" "<<vel(i-1,j-1,0,1)<<"\n";
	    	//   Print()<<"vel(-1,0) "<<vel(i-1,j,0,0)<<" "<<vel(i-1,j,0,1)<<"\n";
	    	//   Print()<<"vel(0,0) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n";
	    	//   Print()<<"vel(0,-1) "<<vel(i,j-1,0,0)<<" "<<vel(i,j-1,0,1)<<"\n\n";
	    	// }
	    	// if (i==0 && j==47){
	    	//   Print()<<" i,j: "<<i<< ","<<j<<"\n";
	    	//   Print()<<"vel(-1,0) "<<vel(i-1,j,0,0)<<" "<<vel(i-1,j,0,1)<<"\n";
	    	//   Print()<<"vel(-1,+1) "<<vel(i-1,j+1,0,0)<<" "<<vel(i-1,j+1,0,1)<<"\n";
	    	//   Print()<<"vel(0,+1) "<<vel(i,j+1,0,0)<<" "<<vel(i,j+1,0,1)<<"\n";
	    	//   Print()<<"vel(0,0) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n";
	    	// }
	    	// if (i==15 && j==47){
	    	//   Print()<<" i,j: "<<i<< ","<<j<<"\n";
	    	//   Print()<<"vel(0,0) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n";
	    	//   Print()<<"vel(0,+1) "<<vel(i,j+1,0,0)<<" "<<vel(i,j+1,0,1)<<"\n";
	    	//   Print()<<"vel(+1,+1) "<<vel(i+1,j+1,0,0)<<" "<<vel(i+1,j+1,0,1)<<"\n";
	    	//   Print()<<"vel(+1,0) "<<vel(i+1,j,0,0)<<" "<<vel(i+1,j,0,1)<<"\n\n";
	    	// }
	     //  }
	    // IntVect mypt={1,35};
	    // if (box.contains(mypt)) {
	    // 	if (i==0 && j==32){
	    // 	  Print()<<"vel(0-1,32-1) "<<vel(i-1,j-1,0,0)<<" "<<vel(i-1,j-1,0,1)<<"\n";
	    // 	  Print()<<"vel(0,32) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n\n";
	    // 	}
	    // 	if (i==0 && j==47){
	    // 	  Print()<<"vel(0-1,47+1) "<<vel(i-1,j+1,0,0)<<" "<<vel(i-1,j+1,0,1)<<"\n";
	    // 	  Print()<<"vel(0,47) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n\n";
	    // 	}
	    //   }
	    // if (i==16 && j==47){
	    //   Print()<<"vel(16-1,47+1) "<<vel(i-1,j+1,0,0)<<" "<<vel(i-1,j+1,0,1)<<"\n";
	    //   Print()<<"vel(0,47) "<<vel(i,j,0,0)<<" "<<vel(i,j,0,1)<<"\n\n";
	    // }

        }
    }
}

AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void mltensor_cross_terms (Box const& box, Array4<Real> const& Ax,
                           Array4<Real const> const& fx,
                           Array4<Real const> const& fy,
                           GpuArray<Real,AMREX_SPACEDIM> const& dxinv,
                           Real bscalar) noexcept
{
    const Real dxi = bscalar * dxinv[0];
    const Real dyi = bscalar * dxinv[1];
    const auto lo = amrex::lbound(box);
    const auto hi = amrex::ubound(box);

    for     (int j = lo.y; j <= hi.y; ++j) {
        AMREX_PRAGMA_SIMD
        for (int i = lo.x; i <= hi.x; ++i) {
            Ax(i,j,0,0) += dxi*(fx(i+1,j  ,0,0) - fx(i,j,0,0))
                +          dyi*(fy(i  ,j+1,0,0) - fy(i,j,0,0));
            Ax(i,j,0,1) += dxi*(fx(i+1,j  ,0,1) - fx(i,j,0,1))
                +          dyi*(fy(i  ,j+1,0,1) - fy(i,j,0,1));
        }
    }
}

}

#endif
