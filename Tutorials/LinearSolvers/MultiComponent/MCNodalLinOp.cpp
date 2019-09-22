//
// Tutorial:    MultiComponent Linear Solve
//
// File:        MCNodalLinOp.cpp
//
// Author:      Brandon Runnels
//              University of Colorado Colorado Springs
//              brunnels@uccs.edu
//              solids.uccs.edu
//
// Date:        September 3, 2019
// 
// Description: This file implements the multi-component linear operator
//
// See also:    ./MCNodalLinOp.H (for definitions and documentation of each function)
//              ./main.cpp (for execution of the tutorial)
// 

#include <AMReX_MLNodeLap_K.H>
#include "MCNodalLinOp.H"

using namespace amrex;

void MCNodalLinOp::Fapply (int amrlev, int mglev, MultiFab& a_out,const MultiFab& a_in) const 
{
	BL_PROFILE("MCNodalLinOp::Fapply()");

	a_out.setVal(0.0);
	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());
	domain.grow(-1); // Shrink domain so we don't operate on any boundaries
	const Real* DX = m_geom[amrlev][mglev].CellSize();

	for (MFIter mfi(a_out, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.tilebox();
		bx.grow(1);        // Expand to cover first layer of ghost nodes
		bx = bx & domain;  // Take intersection of box and the problem domain
			
		amrex::Array4<const amrex::Real> const& in  = a_in.array(mfi);
		amrex::Array4<amrex::Real>       const& out = a_out.array(mfi);

		for (int n = 0; n < getNComp(); n++)
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				out(i,j,k,n) = 0.0;
				for (int m = 0; m < getNComp(); m++)
				{ 
					out(i,j,k,n) += coeff[ncomp*n + m] *
					( AMREX_D_TERM(- (in(i-1,j,k,m) - 2.0 * in(i,j,k,m) + in(i+1,j,k,m)) / DX[0] / DX[0],
					  			   - (in(i,j-1,k,m) - 2.0 * in(i,j,k,m) + in(i,j+1,k,m)) / DX[1] / DX[1],
					  			   - (in(i,j,k-1,m) - 2.0 * in(i,j,k,m) + in(i,j,k+1,m)) / DX[2] / DX[2])) ;
				}

			});
	}
}
void MCNodalLinOp::Diag (int amrlev, int mglev, MultiFab& a_diag)  
{
	a_diag.setVal(1.0);
	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());
	domain.grow(-1); // Shrink domain so we don't operate on any boundaries
	const Real* DX = m_geom[amrlev][mglev].CellSize();

	for (MFIter mfi(a_diag, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.tilebox();
		bx.grow(1);        // Expand to cover first layer of ghost nodes
		bx = bx & domain;  // Take intersection of box and the problem domain
			
		amrex::Array4<amrex::Real> const& diag  = a_diag.array(mfi);

		for (int n = 0; n < getNComp(); n++)
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
				diag(i,j,k,n) = coeff[ncomp*n + n] *
					( AMREX_D_TERM(+ 2.0 / DX[0] / DX[0],
					  			   + 2.0 / DX[1] / DX[1],
					  			   + 2.0 / DX[2] / DX[2]) );
			});
	}
}

void MCNodalLinOp::Diagonal (bool recompute)
{
	if ( !recompute && m_diagonal_computed ) return;
	m_diagonal_computed = true;

	for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	{
		for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		{
			Diag(amrlev,mglev,*m_diag[amrlev][mglev]);
		}
	}
}

void MCNodalLinOp::Fsmooth (int amrlev, int mglev, amrex::MultiFab& a_x, const amrex::MultiFab& a_b) const
{
	BL_PROFILE("MCNodalLinOp::Fsmooth()");

	amrex::Box domain(m_geom[amrlev][mglev].Domain());
	domain.convert(amrex::IntVect::TheNodeVector());
	domain.grow(-1); // Shrink domain so we don't operate on any boundaries

	//int ncomp  = getNComp();
	int nghost = getNGrow(); 

	Real omega = 2./3.; // Damping factor (very important!)

	MultiFab _Ax(a_x.boxArray(), a_x.DistributionMap(), ncomp, nghost);
	MultiFab _Dx(a_x.boxArray(), a_x.DistributionMap(), ncomp, nghost);
	MultiFab _Rx(a_x.boxArray(), a_x.DistributionMap(), ncomp, nghost);
	
	if (!m_diagonal_computed) amrex::Abort("MCNodalLinOp::Diagonal() must be called before using Fsmooth");

	// This is a JACOBI iteration, not Gauss-Seidel.
	// So we need to do twice the number of iterations to get the same behavior as GS.
	for (int ctr = 0; ctr < 2; ctr++)
	{
		Fapply(amrlev,mglev,_Ax,a_x); // find Ax

		amrex::MultiFab::Copy(_Dx,a_x,0,0,ncomp,nghost); // Dx = x
		amrex::MultiFab::Multiply(_Dx,*m_diag[amrlev][mglev],0,0,ncomp,nghost); // Dx *= diag  (Dx = x*diag)

		amrex::MultiFab::Copy(_Rx,_Ax,0,0,ncomp,nghost);     // Rx = Ax
		amrex::MultiFab::Subtract(_Rx,_Dx,0,0,ncomp,nghost); // Rx -= Dx  (Rx = Ax - Dx)


		//for (MFIter mfi(a_x, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
		for (MFIter mfi(a_x, false); mfi.isValid(); ++mfi)
		{
			//Box bx = mfi.tilebox();
			Box bx = mfi.validbox();
			bx.grow(1);        // Expand to cover first layer of ghost nodes
			bx = bx & domain;  // Take intersection of box and the problem domain
			
			amrex::Array4<amrex::Real>       const& x  = a_x.array(mfi);

			amrex::Array4<const amrex::Real> const& b  = a_b.array(mfi);
			amrex::Array4<const amrex::Real> const& Rx =  _Rx.array(mfi);
			amrex::Array4<const amrex::Real> const& diag =  (*m_diag[amrlev][mglev]).array(mfi);

			for (int n = 0; n < getNComp(); n++)
				amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					x(i,j,k,n) = (1.-omega)*x(i,j,k,n) + omega*(b(i,j,k,n) - Rx(i,j,k,n))/diag(i,j,k,n);					
				});
		}
	}
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(a_x,geom);
	nodalSync(amrlev, mglev, a_x);
}

void MCNodalLinOp::normalize (int amrlev, int mglev, MultiFab& a_x) const
{
	BL_PROFILE("MCNodalLinOp::normalize()");
	int nghost = 1; 
	amrex::MultiFab::Divide(a_x,*m_diag[amrlev][mglev],0,0,ncomp,nghost); // Dx *= diag  (Dx = x*diag)
}

 MCNodalLinOp::~MCNodalLinOp () {}

void MCNodalLinOp::define (const Vector<Geometry>& a_geom,
		       const Vector<BoxArray>& a_grids,
		       const Vector<DistributionMapping>& a_dmap,
		       const LPInfo& a_info,
		       const Vector<FabFactory<FArrayBox> const*>& a_factory)
{
	BL_PROFILE("MCNodalLinOp::~Operator()");

	 // This makes sure grids are node-centered;
	 Vector<BoxArray> cc_grids = a_grids;
	 for (auto& ba : cc_grids) {
		 ba.enclosedCells();
	 }

	 MLNodeLinOp::define(a_geom, a_grids, a_dmap, a_info, a_factory);

	 int nghost = 2;
	 // Resize the multifab containing the operator diagonal
	 m_diag.resize(m_num_amr_levels);
	 for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
	 {
		 m_diag[amrlev].resize(m_num_mg_levels[amrlev]);

		 for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
		 {
		 	 m_diag[amrlev][mglev].reset(new MultiFab(amrex::convert(m_grids[amrlev][mglev], amrex::IntVect::TheNodeVector()),
		 						  m_dmap[amrlev][mglev], getNComp(), nghost));
		 }
	 }

	m_lobc.resize(getNComp(),{{AMREX_D_DECL(BCType::Dirichlet,BCType::Dirichlet,BCType::Dirichlet)}});
	m_hibc.resize(getNComp(),{{AMREX_D_DECL(BCType::Dirichlet,BCType::Dirichlet,BCType::Dirichlet)}});
}


void MCNodalLinOp::buildMasks ()
{
	BL_PROFILE("MCNodalLinOp::buildMasks()");
	if (m_masks_built) return;

	m_masks_built = true;

	m_is_bottom_singular = false;

	{
		std::vector< std::pair<int,Box> > isects;
		IArrayBox ccfab;

		for (int amrlev = 0; amrlev < m_num_amr_levels; ++amrlev)
		{
			for (int mglev = 0; mglev < m_num_mg_levels[amrlev]; ++mglev)
			{
				const Geometry& geom = m_geom[amrlev][mglev];
				const auto& period = geom.periodicity();
				const Box& ccdomain = geom.Domain();
				const std::vector<IntVect>& pshifts = period.shiftIntVect();

				Box ccdomain_p = ccdomain;
				for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
					if (geom.isPeriodic(idim)) {
						ccdomain_p.grow(idim, 1);
					}
				}

				{
					auto& dmask = *m_dirichlet_mask[amrlev][mglev];
					const BoxArray& ccba = m_grids[amrlev][mglev];

					for (MFIter mfi(dmask, MFItInfo().SetDynamic(true)); mfi.isValid(); ++mfi)
					{
						const Box& ndbx = mfi.validbox();
						const Box& ccbx = amrex::enclosedCells(ndbx);
						const Box& ccbxg1 = amrex::grow(ccbx,1);
                        
						ccfab.resize(ccbxg1);
						ccfab.setVal(1);
						ccfab.setComplement(2,ccdomain_p,0,1);

						for (const auto& iv : pshifts)
						{
							ccba.intersections(ccbxg1+iv, isects);
							for (const auto& is : isects)
							{
								ccfab.setVal(0, is.second-iv, 0, 1);
							}
						}
					}
				}
			}
		}
	}

	for (int amrlev = 0; amrlev < m_num_amr_levels-1; ++amrlev)
	{
		iMultiFab& cc_mask = *m_cc_fine_mask[amrlev];
		iMultiFab& nd_mask = *m_nd_fine_mask[amrlev];
		LayoutData<int>& has_cf = *m_has_fine_bndry[amrlev];
		const BoxArray& fba = m_grids[amrlev+1][0];
		const BoxArray& cfba = amrex::coarsen(fba, AMRRefRatio(amrlev));

		const Box& ccdom = m_geom[amrlev][0].Domain();

		AMREX_ALWAYS_ASSERT_WITH_MESSAGE(AMRRefRatio(amrlev) == 2, "ref_ratio != 0 not supported");

		cc_mask.setVal(0);  // coarse by default

		const std::vector<IntVect>& pshifts = m_geom[amrlev][0].periodicity().shiftIntVect();

		{
			std::vector< std::pair<int,Box> > isects;

			for (MFIter mfi(cc_mask); mfi.isValid(); ++mfi)
			{
				has_cf[mfi] = 0;
				IArrayBox& fab = cc_mask[mfi];
				const Box& bx = fab.box();
				for (const auto& iv : pshifts)
				{
					cfba.intersections(bx+iv, isects);
					for (const auto& is : isects)
					{
						fab.setVal(1, is.second-iv, 0, 1);
					}
					if (!isects.empty()) has_cf[mfi] = 1;
				}

				amrex_mlndlap_fillbc_cc_i(BL_TO_FORTRAN_ANYD(fab),
							  BL_TO_FORTRAN_BOX(ccdom),
							  m_lobc.data(), m_hibc.data());
			}
		}

		for (MFIter mfi(nd_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
                        Array4<int> const& nmsk = nd_mask.array(mfi);
                        Array4<int const> const& cmsk = cc_mask.const_array(mfi);
                        amrex::ParallelFor(bx, [=] (int i, int j, int k) noexcept
                        {
                            mlndlap_set_nodal_mask(i,j,k,nmsk,cmsk);
                        });
		}
	}

	auto& has_cf = *m_has_fine_bndry[m_num_amr_levels-1];

	for (MFIter mfi(has_cf); mfi.isValid(); ++mfi)
	{
		has_cf[mfi] = 0;
	}

	{
		int amrlev = 0;
		int mglev = m_num_mg_levels[amrlev]-1;
		const iMultiFab& omask = *m_owner_mask[amrlev][mglev];
		m_bottom_dot_mask.define(omask.boxArray(), omask.DistributionMap(), 1, 0);

		const Geometry& geom = m_geom[amrlev][mglev];
		Box nddomain = amrex::surroundingNodes(geom.Domain());

		for (MFIter mfi(m_bottom_dot_mask,true); mfi.isValid(); ++mfi)
		{
			const Box& bx = mfi.tilebox();
                        Array4<Real> const& dfab = m_bottom_dot_mask.array(mfi);
                        Array4<int const> const& sfab = omask.const_array(mfi);
                        mlndlap_set_dot_mask(bx, dfab, sfab, nddomain, m_lobc[0], m_hibc[0]);
		}
	}
}

void MCNodalLinOp::fixUpResidualMask (int amrlev, iMultiFab& resmsk)
{
	BL_PROFILE("MCNodalLinOp::fixUpResidualMask()");

	if (!m_masks_built) buildMasks();

	const iMultiFab& cfmask = *m_nd_fine_mask[amrlev];

#ifdef _OPENMP
#pragma omp parallel
#endif
	for (MFIter mfi(resmsk,true); mfi.isValid(); ++mfi)
	{
            const Box& bx = mfi.tilebox();
            Array4<int> const& rmsk = resmsk.array(mfi);
            Array4<int const> const& fmsk = cfmask.const_array(mfi);
            amrex::ParallelFor(bx, [=] (int i, int j, int k) noexcept
            {
                if (fmsk(i,j,k) == crse_fine_node) rmsk(i,j,k) = 1;
            });
	}

}

void MCNodalLinOp::prepareForSolve ()
{
	BL_PROFILE("MCNodalLinOp::prepareForSolve()");
	MLNodeLinOp::prepareForSolve();
	buildMasks();
	Diagonal(true);
}

void MCNodalLinOp::restriction (int amrlev, int cmglev, MultiFab& crse, MultiFab& fine) const
{
	BL_PROFILE("MCNodalLinOp::restriction()");

	applyBC(amrlev, cmglev-1, fine, BCMode::Homogeneous, StateMode::Solution);

	amrex::Box cdomain = m_geom[amrlev][cmglev].Domain();
	cdomain.convert(amrex::IntVect::TheNodeVector());
	cdomain.grow(-1);

	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), fine.nComp(), fine.nGrow());
	}

	MultiFab* pcrse = (need_parallel_copy) ? &cfine : &crse;

	for (MFIter mfi(*pcrse, amrex::TilingIfNotGPU()); mfi.isValid(); ++mfi)
	{
		Box bx = mfi.tilebox();
		bx = bx & cdomain;

		pcrse->setVal(0.0);
		amrex::Array4<const amrex::Real> const& fdata = fine.array(mfi);
		amrex::Array4<amrex::Real> const& cdata       = pcrse->array(mfi);

		for (int n = 0; n < crse.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=2*I, j=2*J, k=2*K;
						cdata(I,J,K,n) =
							(fdata(i-1,j-1,k-1,n) + fdata(i-1,j-1,k+1,n) + fdata(i-1,j+1,k-1,n) + fdata(i-1,j+1,k+1,n) +
							 fdata(i+1,j-1,k-1,n) + fdata(i+1,j-1,k+1,n) + fdata(i+1,j+1,k-1,n) + fdata(i+1,j+1,k+1,n)) / 64.0
							+
							(fdata(i,j-1,k-1,n) + fdata(i,j-1,k+1,n) + fdata(i,j+1,k-1,n) + fdata(i,j+1,k+1,n) +
							 fdata(i-1,j,k-1,n) + fdata(i+1,j,k-1,n) + fdata(i-1,j,k+1,n) + fdata(i+1,j,k+1,n) +
							 fdata(i-1,j-1,k,n) + fdata(i-1,j+1,k,n) + fdata(i+1,j-1,k,n) + fdata(i+1,j+1,k,n)) / 32.0
							+
							(fdata(i-1,j,k,n) + fdata(i,j-1,k,n) + fdata(i,j,k-1,n) +
							 fdata(i+1,j,k,n) + fdata(i,j+1,k,n) + fdata(i,j,k+1,n)) / 16.0
							+
							fdata(i,j,k,n) / 8.0;
				});
		}
	}

	if (need_parallel_copy) {
		crse.ParallelCopy(cfine);
	}

	amrex::Geometry geom = m_geom[amrlev][cmglev];
	realFillBoundary(crse,geom);
	nodalSync(amrlev, cmglev, crse);
}

void MCNodalLinOp::interpolation (int amrlev, int fmglev, MultiFab& fine, const MultiFab& crse) const
{
	BL_PROFILE("MCNodalLinOp::interpolation()");
	amrex::Box fdomain = m_geom[amrlev][fmglev].Domain(); fdomain.convert(amrex::IntVect::TheNodeVector());
	
	bool need_parallel_copy = !amrex::isMFIterSafe(crse, fine);
	MultiFab cfine;
	const MultiFab* cmf = &crse;
	if (need_parallel_copy) {
		const BoxArray& ba = amrex::coarsen(fine.boxArray(), 2);
		cfine.define(ba, fine.DistributionMap(), crse.nComp(), crse.nGrow());
		cfine.ParallelCopy(crse);
		cmf = &cfine;
	}

	for (MFIter mfi(fine, false); mfi.isValid(); ++mfi)
	{
		const Box& fine_bx = mfi.validbox() & fdomain;
		const Box& course_bx = amrex::coarsen(fine_bx,2);
		const Box& tmpbx = amrex::refine(course_bx,2);
		FArrayBox tmpfab;
		tmpfab.resize(tmpbx,fine.nComp());
		tmpfab.setVal(0.0);
		const amrex::FArrayBox &crsefab = (*cmf)[mfi];

		amrex::Array4<const amrex::Real> const& cdata = crsefab.const_array();
		amrex::Array4<amrex::Real> const& fdata       = tmpfab.array();

		for (int n = 0; n < crse.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (fine_bx,[=] AMREX_GPU_DEVICE(int i, int j, int k) {
					
					int I=i/2, J=j/2, K=k/2;

					if (i%2 == 0 && j%2 == 0 && k%2 ==0) // Coincident
						fdata(i,j,k,n) = cdata(I,J,K,n);
					else if (j%2 == 0 && k%2 == 0) // X Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I+1,J,K,n));
					else if (k%2 == 0 && i%2 == 0) // Y Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I,J+1,K,n));
					else if (i%2 == 0 && j%2 == 0) // Z Edge
						fdata(i,j,k,n) = 0.5 * (cdata(I,J,K,n) + cdata(I,J,K+1,n)); 
					else if (i%2 == 0) // X Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I,J+1,K,n) +
									 cdata(I,J,K+1,n) + cdata(I,J+1,K+1,n));
					else if (j%2 == 0) // Y Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I,J,K+1,n) +
									 cdata(I+1,J,K,n) + cdata(I+1,J,K+1,n));
					else if (k%2 == 0) // Z Face
						fdata(i,j,k,n) = 0.25 * (cdata(I,J,K,n)   + cdata(I+1,J,K,n) +
									 cdata(I,J+1,K,n) + cdata(I+1,J+1,K,n));
					else // Center
						fdata(i,j,k,n) = 0.125 * (cdata(I,J,K,n) +
									  cdata(I+1,J,K,n)   + cdata(I,J+1,K,n)   + cdata(I,J,K+1,n) +
									  cdata(I,J+1,K+1,n) + cdata(I+1,J,K+1,n) + cdata(I+1,J+1,K,n) +
									  cdata(I+1,J+1,K+1,n));
				});
		}
		fine[mfi].plus(tmpfab,fine_bx,fine_bx,0,0,fine.nComp());
	}
	amrex::Geometry geom = m_geom[amrlev][fmglev];
	realFillBoundary(fine,geom);
	nodalSync(amrlev, fmglev, fine);
}
  
void MCNodalLinOp::averageDownSolutionRHS (int camrlev, MultiFab& crse_sol, MultiFab&,
				  const MultiFab& fine_sol, const MultiFab&)
{
	BL_PROFILE("MCNodalLinOp::averageDownSolutionRHS()");
	const auto& amrrr = AMRRefRatio(camrlev);
	amrex::average_down(fine_sol, crse_sol, 0, crse_sol.nComp(), amrrr);
	if (isSingular(0)) amrex::Abort("Singular operators not supported!");
}

void MCNodalLinOp::realFillBoundary(MultiFab &phi, const Geometry &geom) const
{
	for (int i = 0; i < 2; i++)
	{
		MultiFab & mf = phi;
		mf.FillBoundary(geom.periodicity());
		//const int ncomp = mf.nComp();
		const int ng1 = 1;
		const int ng2 = 2;
		MultiFab tmpmf(mf.boxArray(), mf.DistributionMap(), ncomp, ng1);
		MultiFab::Copy(tmpmf, mf, 0, 0, ncomp, ng1); 
		mf.ParallelCopy   (tmpmf, 0, 0, ncomp, ng1, ng2, geom.periodicity());
	}
}

void MCNodalLinOp::applyBC (int amrlev, int mglev, MultiFab& phi, BCMode,
		   					amrex::MLLinOp::StateMode , bool skip_fillboundary) const
{
	BL_PROFILE("MCNodalLinOp::applyBC()");
	const Geometry& geom = m_geom[amrlev][mglev];
	if (!skip_fillboundary) realFillBoundary(phi,geom);
}

void MCNodalLinOp::reflux (int crse_amrlev,
		  MultiFab& res, const MultiFab&, const MultiFab&,
		  MultiFab& fine_res, MultiFab&, const MultiFab&) const
{
	BL_PROFILE("MCNodalLinOp::reflux()");

	amrex::Box cdomain(m_geom[crse_amrlev][0].Domain());
	cdomain.convert(amrex::IntVect::TheNodeVector());

	const Geometry& cgeom = m_geom[crse_amrlev  ][0];

 	const BoxArray&            fba = fine_res.boxArray();
 	const DistributionMapping& fdm = fine_res.DistributionMap();

 	MultiFab fine_res_for_coarse(amrex::coarsen(fba, 2), fdm, ncomp, 2);
	fine_res_for_coarse.ParallelCopy(res,0,0,ncomp,0,0,cgeom.periodicity());

 	applyBC(crse_amrlev+1, 0, fine_res, BCMode::Inhomogeneous, StateMode::Solution);

	const int coarse_fine_node = 1;
	const int fine_fine_node = 2;

	amrex::iMultiFab nodemask(amrex::coarsen(fba,2), fdm, 1, 2);
	nodemask.ParallelCopy(*m_nd_fine_mask[crse_amrlev],0,0,1,0,0,cgeom.periodicity());

	amrex::iMultiFab cellmask(amrex::convert(amrex::coarsen(fba,2),amrex::IntVect::TheCellVector()), fdm, 1, 2);
	cellmask.ParallelCopy(*m_cc_fine_mask[crse_amrlev],0,0,1,1,1,cgeom.periodicity());
	
	for (MFIter mfi(fine_res_for_coarse, false); mfi.isValid(); ++mfi)
	{
		const Box& bx = mfi.validbox();

		amrex::Array4<const int> const& nmask = nodemask.array(mfi);
		amrex::Array4<amrex::Real> const& cdata = fine_res_for_coarse.array(mfi);
		amrex::Array4<const amrex::Real> const& fdata       = fine_res.array(mfi);

		const Dim3 lo= amrex::lbound(cdomain), hi = amrex::ubound(cdomain);

		for (int n = 0; n < fine_res.nComp(); n++)
		{
			// I,J,K == coarse coordinates
			// i,j,k == fine coordinates
			amrex::ParallelFor (bx,[=] AMREX_GPU_DEVICE(int I, int J, int K) {
					int i=I*2, j=J*2, k=K*2;
					
					if (nmask(I,J,K) == fine_fine_node || nmask(I,J,K) == coarse_fine_node)
						{
							if ((I == lo.x || I == hi.x) &&
							    (J == lo.y || J == hi.y) &&
							    (K == lo.z || K == hi.z)) // Corner
								cdata(I,J,K,n) = fdata(i,j,k,n);
							else if ((J == lo.y || J == hi.y) &&
								 (K == lo.z || K == hi.z)) // X edge
								cdata(I,J,K,n) = 0.25*fdata(i-1,j,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i+1,j,k,n);
							else if ((K == lo.z || K == hi.z) &&
								 (I == lo.x || I == hi.x)) // Y edge
								cdata(I,J,K,n) = 0.25*fdata(i,j-1,k,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j+1,k,n);
							else if ((I == lo.x || I == hi.x) &&
								 (J == lo.y || J == hi.y)) // Z edge
								cdata(I,J,K,n) = 0.25*fdata(i,j,k-1,n) + 0.5*fdata(i,j,k,n) + 0.25*fdata(i,j,k+1,n);
							else if (I == lo.x || I == hi.x) // X face
								cdata(I,J,K,n) =
									(+     fdata(i,j-1,k-1,n) + 2.0*fdata(i,j,k-1,n) +     fdata(i,j+1,k-1,n)
									 + 2.0*fdata(i,j-1,k  ,n) + 4.0*fdata(i,j,k  ,n) + 2.0*fdata(i,j+1,k  ,n) 
									 +     fdata(i,j-1,k+1,n) + 2.0*fdata(i,j,k+1,n) +     fdata(i,j+1,k+1,n))/16.0;
							else if (J == lo.y || J == hi.y) // Y face
								cdata(I,J,K,n) =
									(+     fdata(i-1,j,k-1,n) + 2.0*fdata(i-1,j,k,n) +     fdata(i-1,j,k+1,n)
									 + 2.0*fdata(i  ,j,k-1,n) + 4.0*fdata(i  ,j,k,n) + 2.0*fdata(i  ,j,k+1,n) 
									 +     fdata(i+1,j,k-1,n) + 2.0*fdata(i+1,j,k,n) +     fdata(i+1,j,k+1,n))/16.0;
							else if (K == lo.z || K == hi.z) // Z face
								cdata(I,J,K,n) =
									(+     fdata(i-1,j-1,k,n) + 2.0*fdata(i,j-1,k,n) +     fdata(i+1,j-1,k,n)
									 + 2.0*fdata(i-1,j  ,k,n) + 4.0*fdata(i,j  ,k,n) + 2.0*fdata(i+1,j  ,k,n) 
									 +     fdata(i-1,j+1,k,n) + 2.0*fdata(i,j+1,k,n) +     fdata(i+1,j+1,k,n))/16.0;
							else // Interior
								cdata(I,J,K,n) =
									(fdata(i-1,j-1,k-1,n) + fdata(i-1,j-1,k+1,n) + fdata(i-1,j+1,k-1,n) + fdata(i-1,j+1,k+1,n) +
									 fdata(i+1,j-1,k-1,n) + fdata(i+1,j-1,k+1,n) + fdata(i+1,j+1,k-1,n) + fdata(i+1,j+1,k+1,n)) / 64.0
									+
									(fdata(i,j-1,k-1,n) + fdata(i,j-1,k+1,n) + fdata(i,j+1,k-1,n) + fdata(i,j+1,k+1,n) +
									 fdata(i-1,j,k-1,n) + fdata(i+1,j,k-1,n) + fdata(i-1,j,k+1,n) + fdata(i+1,j,k+1,n) +
									 fdata(i-1,j-1,k,n) + fdata(i-1,j+1,k,n) + fdata(i+1,j-1,k,n) + fdata(i+1,j+1,k,n)) / 32.0
									+
									(fdata(i-1,j,k,n) + fdata(i,j-1,k,n) + fdata(i,j,k-1,n) +
									 fdata(i+1,j,k,n) + fdata(i,j+1,k,n) + fdata(i,j,k+1,n)) / 16.0
									+
									fdata(i,j,k,n) / 8.0;
						}

				});
		}
	}

	// Copy the fine residual restricted onto the coarse grid
	// into the final residual.
	res.ParallelCopy(fine_res_for_coarse,0,0,ncomp,0,0,cgeom.periodicity());

	const int mglev = 0;

	// Sync up ghost nodes
	amrex::Geometry geom = m_geom[crse_amrlev][mglev];
	realFillBoundary(res,geom);
	nodalSync(crse_amrlev,mglev, res);
	return;
}

void
MCNodalLinOp::solutionResidual (int amrlev, MultiFab& resid, MultiFab& x, const MultiFab& b,
			    const MultiFab*)
{
	const int mglev = 0;
	apply(amrlev, mglev, resid, x, BCMode::Inhomogeneous, StateMode::Solution);
	MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, 2);
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(resid,geom);
}

void
MCNodalLinOp::correctionResidual (int amrlev, int mglev, MultiFab& resid, MultiFab& x, const MultiFab& b,
			      BCMode, const MultiFab* )
{
	resid.setVal(0.0);
	apply(amrlev, mglev, resid, x, BCMode::Homogeneous, StateMode::Correction);
	MultiFab::Xpay(resid, -1.0, b, 0, 0, ncomp, resid.nGrow());
	amrex::Geometry geom = m_geom[amrlev][mglev];
	realFillBoundary(resid,geom);
}
