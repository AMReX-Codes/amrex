
#include <REAL.H>

extern "C" 
{
    void pxr_set_gamma(const long* np,
		       const Real* up, const Real* vp, const Real* wp,
		       const Real* gaminv);

    // void depose_rho_vecHVv2_1_1_1( Real* rho, const long* np,
    // 				   const Real* xp, const Real* yp, const Real* zp,
    // 				   const Real* w, const Real* q,
    // 				   const Real* xmin, const Real* ymin, const Real* zmin,
    // 				   const Real* dx, const Real* dy, const Real* dz,
    // 				   const long* nx, const long* ny, const long* nz,
    // 				   const long* nxguard, const long* nyguard, const long* nzguard,
    // 				   const long* lvect );
    
    // void depose_rho_scalar_1_1_1( Real* rho, const long* np,
    // 				  const Real* xp, const Real* yp, const Real* zp,
    // 				  const Real* w, const Real* q,
    // 				  const Real* xmin, const Real* ymin, const Real* zmin,
    // 				  const Real* dx, const Real* dy, const Real* dz,
    // 				  const long* nx, const long* ny, const long* nz,
    // 				  const long* nxguard, const long* nyguard, const long* nzguard,
    // 				  const long* lvect );

    // void depose_rho_scalar_2_2_2( Real* rho, const long* np,
    // 				  const Real* xp, const Real* yp, const Real* zp,
    // 				  const Real* w, const Real* q,
    // 				  const Real* xmin, const Real* ymin, const Real* zmin,
    // 				  const Real* dx, const Real* dy, const Real* dz,
    // 				  const long* nx, const long* ny, const long* nz,
    // 				  const long* nxguard, const long* nyguard, const long* nzguard,
    // 				  const long* lvect );

    // void depose_rho_scalar_3_3_3( Real* rho, const long* np,
    // 				  const Real* xp, const Real* yp, const Real* zp,
    // 				  const Real* w, const Real* q,
    // 				  const Real* xmin, const Real* ymin, const Real* zmin,
    // 				  const Real* dx, const Real* dy, const Real* dz,
    // 				  const long* nx, const long* ny, const long* nz,
    // 				  const long* nxguard, const long* nyguard, const long* nzguard,
    // 				  const long* lvect );


    void depose_jxjyjz_vecHVv2_1_1_1(Real* jx, Real* jy, Real* jz, const long* np,
				     const Real* xp, const Real* yp, const Real* zp,
				     const Real* uxp, const Real* uyp,const Real* uzp,
				     const Real* gip, const Real* w, const Real* q,
				     const Real* xmin, const Real* ymin, const Real* zmin,
				     const Real* dt, 
				     const Real* dx, const Real* dy, const Real* dz,
				     const long* nx, const long* ny, const long* nz,
				     const long* nxguard, const long* nyguard, const long* nzguard);

    void depose_jxjyjz_scalar_1_1_1(Real* jx, Real* jy, Real* jz, const long* np,
				    const Real* xp, const Real* yp, const Real* zp,
				    const Real* uxp, const Real* uyp,const Real* uzp,
				    const Real* gip, const Real* w, const Real* q,
				    const Real* xmin, const Real* ymin, const Real* zmin,
				    const Real* dt, 
				    const Real* dx, const Real* dy, const Real* dz,
				    const long* nx, const long* ny, const long* nz,
				    const long* nxguard, const long* nyguard, const long* nzguard);

    void warpx_geteb3d_energy_conserving(const long* np, 
					 const Real* xp, const Real* yp, const Real* zp,
					 const Real* exp, const Real* eyp,const Real* ezp,
					 const Real* bxp, const Real* byp,const Real* bzp,
					 const Real* xmin, const Real* ymin, const Real* zmin,
					 const Real* dx, const Real* dy, const Real* dz,
					 const long* nx, const long* ny, const long* nz,
					 const long* nxguard, const long* nyguard, const long* nzguard,
					 const long* nox, const long* noy, const long* noz,
					 const Real* exg, const Real* eyg, const Real* ezg,
					 const Real* bxg, const Real* byg, const Real* bzg,
					 const int* ll4symtry, const int* l_lower_order_in_v,  
					 const long* field_gathe_algo);

    void pxr_epush_v(const long* np, 
		     const Real* uxp, const Real* uyp, const Real* uzp,
		     const Real* exp, const Real* eyp,const Real* ezp,
		     const Real* charge, const Real* mass, Real* dt);
    
    void pxr_bpush_v(const long* np, 
		     const Real* uxp, const Real* uyp, const Real* uzp,
		     const Real* gaminv, 
		     const Real* bxp, const Real* byp,const Real* bzp,
		     const Real* charge, const Real* mass, Real* dt);
    
    void pxr_pushxyz(const long* np, 
		     const Real*  xp, const Real*  yp, const Real*  zp,
		     const Real* uxp, const Real* uyp, const Real* uzp,
		     const Real* gaminv, Real* dt);

    void warpx_pxrpush_em3d_evec_norder( Real* ex, Real* ey, Real* ez,
					 Real* bx, Real* by, Real* bz,
					 const Real* jx, const Real* jy, const Real* jz,
					 const Real* mudt,
					 const Real* dtsdx, const Real* dtsdy, const Real* dtsdz, 
					 const long* nxlocal, const long* nylocal, const long* nzlocal,
					 const long* nxorder, const long* nyorder, const long* nzorder,
					 const long* nxguard, const long* nyguard, const long* nzguard, 
					 const long* nxs, const long* nys, const long* nzs,
					 const int* l_nodal );
    
    void warpx_pxrpush_em3d_bvec_norder( Real* ex, Real* ey, Real* ez,
					 Real* bx, Real* by, Real* bz,
					 const Real* dtsdx, const Real* dtsdy, const Real* dtsdz, 
					 const long* nxlocal, const long* nylocal, const long* nzlocal,
					 const long* nxorder, const long* nyorder, const long* nzorder,
					 const long* nxguard, const long* nyguard, const long* nzguard, 
					 const long* nxs, const long* nys, const long* nzs,
					 const int* l_nodal );
}
