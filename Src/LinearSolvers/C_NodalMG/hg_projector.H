//BL_COPYRIGHT_NOTICE

#ifndef _HG_PROJECTOR_H_
#define _HG_PROJECTOR_H_

#include "hg_multi.H"
#include <Geometry.H>

extern PArray<MultiFab> null_amr_real;

class holy_grail_amr_projector : public holy_grail_amr_multigrid 
{
public:
    //
    // No, I'm not pulling your leg.  Lev_min and Lev_max, supplied with
    // project, tell the range of levels over which to project.  The
    // arguments to the constructor give the permissible range for Lev_min
    // and Lev_max.  Lev_min can be anything from Lev_min_min to Lev_min_max;
    // Lev_max can then be anything from Lev_min to Lev_max_max.
    //
    holy_grail_amr_projector (const Array<BoxArray>&          Mesh,
                              const Array<IntVect>&           Gen_ratio,
                              const Box&                      fdomain,
                              int                             Lev_min_min,
                              int                             Lev_min_max,
                              int                             Lev_max_max,
                              const amr_fluid_boundary_class& Boundary,
                              stencil                         stencil_ = cross,
                              int                             Pcode = 2)
	:
        holy_grail_amr_multigrid(Mesh,Gen_ratio,fdomain,Lev_min_min,
                                 Lev_min_max,Lev_max_max,Boundary,
                                 stencil_,
                                 Pcode),
        make_sparse_node_source_solvable(false)
    {}     
    const Array<BoxArray>& mesh () const
    {
	return ml_mesh;
    }
    //
    // All of the following projection calls take a vector field u and
    // a scalar field p.  On return p will contain the computed potential
    // which satisfies Poisson's equation.  The field contained in p on
    // entry to these routines will be used as an initial guess for the
    // iterative scheme.  Also on return u, if it exists, will be updated
    // to u - Gp.
    //
    // For project, the right hand side will be Du.
    // For sync_project, the right hand side will be Du computed at the
    // coarse-fine interface, and 0 everywhere in the grid interiors.
    // For manual_project, the right hand side is a separate input field,
    // and u can be passed in 0(NULL) if you don't care about the update.
    // If use_u is false (0), the right hand side will be used without
    // modification.  If it is true (1), Du will be added to the rhs provided.
    
    // Let start_norm = max[norm(div u), norm(first residual)].
    // Project will attempt to reduce the Holy Grail residual to
    // less than max[tol * start_norm, scale].  At least one multilevel
    // cycle will be executed down to lev_min.  The last cycle will omit
    // any actions on or below lev_min.
    
    // H is an array of SPACEDIM Reals containing the fine mesh spacing
    // in each dimension.  "fine" means lev_max.  If a positive value is
    // not supplied for lev_min, lev_min_max will be used.  If
    // a positive value is not supplied for lev_max, lev_min will be used.
    
    // Ghost cells:  All version require exactly 1 row of ghost cells on p.
    // project and sync_project, which compute a divergence, require
    // exactly 1 row of ghost cells on u (manual_project does not have
    // this requirement).  Sigma does not require ghost cells.  The only
    // case where user-supplied information in the ghost cells is actually
    // used is the normal velocity at an inflow boundary---in all other
    // cases the projection will fill in ghost cells according to its
    // own requirements.
    //
    void project (PArray<MultiFab>* u,
                  PArray<MultiFab>& p,
                  PArray<MultiFab>& Coarse_source,
                  MultiFab*         Sync_resid_crse,
                  MultiFab*         Sync_resid_fine,
                  const Geometry&   crse_geom,
                  Real              H[],
                  Real              tol,
                  int               Lev_min = -1,
                  int               Lev_max = -1,
                  Real              scale = 0.0) 
    {
	project(u, p, Coarse_source, null_amr_real, 
                Sync_resid_crse, Sync_resid_fine, crse_geom, 
                H, tol, Lev_min, Lev_max, scale);
    }
    void project (PArray<MultiFab>* u,
                  PArray<MultiFab>& p,
                  PArray<MultiFab>& Coarse_source,
                  PArray<MultiFab>& Sigma,
                  MultiFab*         Sync_resid_crse,
                  MultiFab*         Sync_resid_fine,
                  const Geometry&   crse_geom,
                  Real              H[],
                  Real              tol,
                  int               Lev_min = -1,
                  int               Lev_max = -1,
                  Real              scale = 0.0);

    void sync_project (PArray<MultiFab>* u,
                       PArray<MultiFab>& p,
                       PArray<MultiFab>& Coarse_source,
                       Real              H[],
                       Real              tol,
                       int               Lev_min = -1,
                       int               Lev_max = -1,
                       Real              scale = 0.0) 
    {
	sync_project(u, p, Coarse_source, null_amr_real, H, tol, Lev_min, Lev_max, scale);
    }

    void sync_project (PArray<MultiFab>* u,
                       PArray<MultiFab>& p,
                       PArray<MultiFab>& Coarse_source,
                       PArray<MultiFab>& Sigma,
                       Real              H[],
                       Real              tol,
                       int               Lev_min = -1,
                       int               Lev_max = -1,
                       Real              scale = 0.0);
    //
    // This is the do-it-yourself version:  The right hand side for the
    // elliptic solve is specified in rhs.  The resulting potential will
    // be returned in p, and its gradient will be subtracted from u
    // if u is non-null.
    //
    // If rhs is node-based it will simply be used as the right hand side.
    // No solvability modifications will be made.
    //
    // If rhs is cell-based it will be conservatively averaged onto nodes
    // to use as the right hand side.  rhs must be passed in with a one
    // cell wide border.  At inflow boundaries border values should
    // be set to correct inflow condition.  Other border values passed
    // in may be meaningless, but should not be NaNs.
    //
    // This routine will modify the borders of the cell-based rhs.  Also,
    // if the problem being solved is singular, rhs will be adjusted so
    // that it integrates to 0 to maximum precision.
    //
    void manual_project (PArray<MultiFab>* u,
                         PArray<MultiFab>& p,
                         PArray<MultiFab>& rhs,
                         PArray<MultiFab>& Coarse_source,
                         MultiFab*         Sync_resid_crse,
                         MultiFab*         Sync_resid_fine,
                         const Geometry&   crse_geom,
                         bool              use_u,
                         Real              H[],
                         Real              tol,
                         int               Lev_min = -1,
                         int               Lev_max = -1,
                         Real              scale = 0.0) 
    {
	manual_project(u, p, rhs, Coarse_source, null_amr_real, 
                       Sync_resid_crse, Sync_resid_fine, crse_geom, 
                       use_u, H, tol, Lev_min, Lev_max, scale);
    }

    void manual_project (PArray<MultiFab>* u,
                         PArray<MultiFab>& p,
                         PArray<MultiFab>& rhs,
                         PArray<MultiFab>& Coarse_source,
                         PArray<MultiFab>& Sigma,
                         MultiFab*         Sync_resid_crse,
                         MultiFab*         Sync_resid_fine,
                         const Geometry&   crse_geom,
                         bool              use_u,
                         Real              H[],
                         Real              tol,
                         int               Lev_min = -1,
                         int               Lev_max = -1,
                         Real              scale = 0.0);

    void fill_sync_reg (PArray<MultiFab>* u_local,
                        PArray<MultiFab>& p,
                        PArray<MultiFab>& rhs_local,
                        PArray<MultiFab>& Sigma_local,
                        MultiFab*         Sync_resid,
                        const Geometry&   crse_geom,
                        Real              H[],
                        int               Lev_min,
                        int               is_fine);
    
    void make_it_so ()
    {
	make_sparse_node_source_solvable = true;
    }

protected:

    void right_hand_side (PArray<MultiFab>* u, PArray<MultiFab>& S,
                          int for_sync_reg);
    void sync_right_hand_side (PArray<MultiFab>* u);
    void interface_average (PArray<MultiFab>& S, int lev);
    void interface_divergence (PArray<MultiFab>* u, int lev);
    void form_solution_vector (PArray<MultiFab>* u, const PArray<MultiFab>& sigma_in);
    void sparse_node_source_adjustment (PArray<MultiFab>& sparse_source);
    void grid_divergence (PArray<MultiFab>* u, PArray<MultiFab>& s,
                          int for_sync_reg);
    void grid_average (PArray<MultiFab>& S, PArray<MultiFab>& src,        
                       int for_sync_reg);
    //
    // The data.
    //
    bool make_sparse_node_source_solvable;    
};

#endif /*_HG_PROJECTOR_H_*/
