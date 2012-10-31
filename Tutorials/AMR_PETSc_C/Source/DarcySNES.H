#ifndef _DarcySNES_H_
#define _DarcySNES_H_

#include <MFTower.H>
#include <BCRec.H>
#include <Layout.H>
#include <MLBoundary.H>

#include <petsc.h>

class DarcySNES
{
public:

  DarcySNES(Layout&     _layout,
            MLBoundary& _mlb);

  ~DarcySNES();

  static void SetVerbose(int v);
  static int GetVerbose();

  void ComputeDarcyFlux(PArray<MFTower>& darcy_flux,
                        int              fComp,
                        MFTower&         pressure,
                        int              pComp,
                        MFTower&         reduced_saturation,
                        int              rsComp);
  
  void ReducedSaturationGivenPressure(const FArrayBox&      pressure,
				      int                   pComp,
				      FArrayBox&            reduced_saturation,
				      int                   satComp,
				      const Box&            box,
				      const Layout::IntFab& mask,
                                      Real                  maskedVal = 0);

  void KrGivenReducedSaturation(const FArrayBox& reduced_saturation,
				int              satComp,
				FArrayBox&       Kr,
				int              kComp,
				const Box&       box,
				const Layout::IntFab& mask,
                                Real                  maskedVal = 0);

  void RhoSatGivenReducedSaturation(const FArrayBox&      reduced_saturation,
                                    int                   redComp,
                                    FArrayBox&            RhoSat,
                                    int                   rsComp,
                                    const Box&            box,
                                    const Layout::IntFab& mask,
                                    Real                  maskedVal = 0);

  void DpDtResidual(MFTower& residual,
		    int      resComp,
		    MFTower& pressure,
		    int      pComp,
		    Real     dt);

  void DivRhoU(MFTower& DivRhoU,
	       int      druComp,
	       MFTower& pressure,
	       int      pComp,
               MFTower& rhoSat,
               int      rsComp);
  
  int Solve(PArray<MultiFab>& RhoSat_new,
	    PArray<MultiFab>& RhoSat_old,
	    int               rhosatComp,
	    PArray<MultiFab>& P_new,
	    int               pComp,
	    Real              dt);
    
  void BuildOpSkel(Mat& J);

  void FillPatch(MFTower& mft,
		 int      dComp,
		 int      nComp,
		 bool     do_piecewise_constant=false);
  
  void SetTheDarcySNES(DarcySNES* ptr);

  Layout& layout;
  MLBoundary& mlb;
  int nLevs;
  MFTFillPatch* mftfp;

  SNES           snes;
  ISColoring     iscoloring;
  MatFDColoring  matfdcoloring;

  // Scratch (large) space
  MFTower *Rhs;
  PArray<MFTower> DarcyFlux;

  Vec RhsV, SolnV, Wrk1V, Wrk2V;
  Mat Jac;

  // Shell structures, pointing to memory owned by someone else
  MFTower *RhoSatOld, *RhoSatNew, *Pnew;
  int pressure_comp, rhoSat_comp; // Corresponding components of this "borrowed" data

  Array<Real> gravity;
  Real density, mu, porosity;
  Real kappax,kappaz,sigma,Sr,m;
  Real dt;

  Array<std::map<int, Array<FArrayBox*> > > bc_pressure_values; // val[lev][gidx][section_i]
  Array< Array<std::map<int, Array<FArrayBox*> > > > bc_flux_values; // flux[d][lev][gidx][section_i]
  Array<Orientation> dirichlet_faces;
  Array<Orientation> neumann_faces;

protected:
  static void CleanupStatics();
  void Initialize();
  void ResetRhoSat(PArray<MultiFab>& S_new,
		   PArray<MultiFab>& S_old,
		   int               rsComp,
		   PArray<MultiFab>& P_new,
		   int               pComp);

  static int verbose;
};

#endif
