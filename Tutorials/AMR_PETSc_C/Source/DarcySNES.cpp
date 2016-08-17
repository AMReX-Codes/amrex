#include <DarcySNES.H>
#include <DarcySNES_F.H>
#include <Utility.H>
#include <Layout.H>
#include <VisMF.H>

#undef PETSC_3_2
#define PETSC_3_2 1

using std::cout;
using std::cerr;
using std::endl;

  // Default material parameters (Soil)
static Real kappax_m2_DEF = 5.13e-13;
static Real kappaz_m2_DEF = 2.87e-13;
static Real sigma_DEF = .000302; 
static Real Sr_DEF = 0.354;   
static Real m_DEF = 0.291;   
static Real se_threshold = 0.9;

static int  pressure_maxorder = 3;
static bool abort_on_nl_fail = false;
static bool centered_diff_J = true;

static Real porosity_DEF = 0.3;
static Real density_kg_m3_DEF = 998.2;
static Real mu_sPa_DEF = .001002;
static Real gravity_mag_DEF = 9.8;
static int  gravity_dir_DEF = BL_SPACEDIM - 1;
static Real Pa_per_atm = 1.01325e5;
static int  verbose_DEF = 0;

static int verbose_SNES = 0;

static bool initialized = false;

int DarcySNES::verbose = verbose_DEF;

void
DarcySNES::CleanupStatics ()
{
  verbose = verbose_DEF;
  initialized = false;
}

void
DarcySNES::Initialize()
{
  BoxLib::ExecOnFinalize(DarcySNES::CleanupStatics);
  initialized = true;
}

void
DarcySNES::SetVerbose(int v)
{
  verbose = v;
}

int
DarcySNES::GetVerbose()
{
  return verbose;
}

bool Vlevel (int verbose, int level) {
  return verbose>level && ParallelDescriptor::IOProcessor();
}

void
Vlevel4 (const std::string& msg, int verbose) {
  if (Vlevel(verbose,4)) {
    std::cout << msg << std::endl;
  }
}

void
Vlevel3 (const std::string& msg, int verbose) {
  if (Vlevel(verbose,3)) {
    std::cout << msg << std::endl;
  }
}

void
Vlevel2 (const std::string& msg, int verbose) {
  if (Vlevel(verbose,2)) {
    std::cout << msg << std::endl;
  }
}

void
Vlevel1 (const std::string& msg, int verbose) {
  if (Vlevel(verbose,1)) {
    std::cout << msg << std::endl;
  }
}

// Forward declaration of local helper functions
PetscErrorCode DarcyComputeJacobianColor(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx);
PetscErrorCode DarcyMatFDColoringApply(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx);
PetscErrorCode DarcyRes_DpDt(SNES snes,Vec x,Vec f,void *dummy);

DarcySNES::DarcySNES(Layout&      _layout,
                     MLBoundary&  _mlb)
  : layout(_layout), mlb(_mlb)
{
  Vlevel1("**** DarcySNES constructor **** (begin)",verbose);

  if (!initialized) {
    Initialize();
  }

  // We will be evaluating the darcy velocity as:
  //
  //  u_Darcy [m/s] = ( kappa [X . mD] / mu [Pa.s] ).Grad(p) [atm/m]
  //
  // where X is the factor necessary to have this formula be dimensionally
  // consistent.  If kappa is in mD, X is 1.e-10.  In fact, kappa is 
  // specified in m^2.  kappa_m2 * 1.01325e5 * 1.e10 = kappa_mD, then we fold
  // the 1.e-10 into this so that u = (kappa/mu) Grad(p) 
  //
  // Similarly, sigma is given in 1/Pa, so we convert to 1/atm
  kappax = kappax_m2_DEF * Pa_per_atm; // kappa in 1.e-10 mD
  kappaz = kappaz_m2_DEF * Pa_per_atm;
  sigma = sigma_DEF * Pa_per_atm; // sigma in 1/atm
  Sr = Sr_DEF;
  m = m_DEF;
  porosity = porosity_DEF;
  density = density_kg_m3_DEF; // rho in kg/m^3
  mu = mu_sPa_DEF; // dynamic viscosity, in Pa.s
  gravity.resize(BL_SPACEDIM,0); gravity[gravity_dir_DEF] = gravity_mag_DEF / Pa_per_atm; // in atm.m^2/kg so rho.g in atm/m

  nLevs = layout.NumLevels();
  mftfp = new MFTFillPatch(layout);

  // These will be set prior to each solve call in order to support the case that
  // the unlying multifabs get changed between repeated calls to this solver (AMR
  // typically advances the state data by swapping the underlying pointers).
  RhoSatOld = 0;
  RhoSatNew = 0;
  Pnew = 0;

  int nComp = 1;
  int nGrow = 1;
  int nGrowFlux = 0;

  Rhs = new MFTower(layout,MFTower::CC,nComp,nGrow,nLevs);
  DarcyFlux.resize(BL_SPACEDIM,PArrayManage);
  for (int d=0; d<BL_SPACEDIM; ++d) {
    DarcyFlux.set(d, new MFTower(layout,MFTower::EC[d],nComp,nGrowFlux,nLevs));
  }

  PetscErrorCode ierr;       
  int n = layout.NumberOfLocalNodeIds();
  int N = layout.NumberOfGlobalNodeIds();
  MPI_Comm comm = ParallelDescriptor::Communicator();
  ierr = VecCreateMPI(comm,n,N,&RhsV); CHKPETSC(ierr);
  ierr = VecDuplicate(RhsV,&SolnV); CHKPETSC(ierr);
  if (centered_diff_J) {
    ierr = VecDuplicate(RhsV,&Wrk1V); CHKPETSC(ierr);
    ierr = VecDuplicate(RhsV,&Wrk2V); CHKPETSC(ierr);
  }

  mftfp->BuildStencil(mlb.BC(), pressure_maxorder);

  // Estmated number of nonzero local columns of J
  int d_nz = 1 + (pressure_maxorder-1)*(2*BL_SPACEDIM); 
  int o_nz = 0; // Estimated number of nonzero nonlocal (off-diagonal) columns of J

#if defined(PETSC_3_2)
  ierr = MatCreateMPIAIJ(comm, n, n, N, N, d_nz, PETSC_NULL, o_nz, PETSC_NULL, &Jac); CHKPETSC(ierr);
#else
  ierr = MatCreate(comm, &Jac); CHKPETSC(ierr);
  ierr = MatSetSizes(Jac,n,n,N,N);  CHKPETSC(ierr);
  ierr = MatSetFromOptions(Jac); CHKPETSC(ierr);
  ierr = MatSeqAIJSetPreallocation(Jac, d_nz*d_nz, PETSC_NULL); CHKPETSC(ierr);
  ierr = MatMPIAIJSetPreallocation(Jac, d_nz, PETSC_NULL, o_nz, PETSC_NULL); CHKPETSC(ierr);
#endif  

  BuildOpSkel(Jac);
  
  matfdcoloring = 0;
  ierr = SNESCreate(comm,&snes); CHKPETSC(ierr);
  ierr = SNESSetFunction(snes,RhsV,DarcyRes_DpDt,(void*)(this)); CHKPETSC(ierr);

  ierr = MatGetColoring(Jac,MATCOLORINGSL,&iscoloring); CHKPETSC(ierr);
  ierr = MatFDColoringCreate(Jac,iscoloring,&matfdcoloring); CHKPETSC(ierr);
  ierr = MatFDColoringSetFunction(matfdcoloring,
                                  (PetscErrorCode (*)(void))DarcyRes_DpDt,
                                  (void*)(this)); CHKPETSC(ierr);
  ierr = MatFDColoringSetFromOptions(matfdcoloring); CHKPETSC(ierr);
  ierr = SNESSetJacobian(snes,Jac,Jac,DarcyComputeJacobianColor,matfdcoloring);CHKPETSC(ierr);
  ierr = SNESSetFromOptions(snes);CHKPETSC(ierr);

  Vlevel1("**** DarcySNES constructor **** (end)",verbose);
}

DarcySNES::~DarcySNES()
{
  Vlevel1("**** DarcySNES destructor **** (begin)",verbose);
  delete Pnew;
  delete RhoSatNew;
  delete RhoSatOld;

  PetscErrorCode ierr;
  ierr = MatFDColoringDestroy(&matfdcoloring); CHKPETSC(ierr);
  ierr = ISColoringDestroy(&iscoloring);
  ierr = SNESDestroy(&snes); CHKPETSC(ierr);
  ierr = MatDestroy(&Jac); CHKPETSC(ierr);
  
  if (centered_diff_J) {
    ierr = VecDestroy(&Wrk1V); CHKPETSC(ierr);
    ierr = VecDestroy(&Wrk2V); CHKPETSC(ierr);
  }
  ierr = VecDestroy(&SolnV); CHKPETSC(ierr);
  ierr = VecDestroy(&RhsV); CHKPETSC(ierr);

  DarcyFlux.clear();  
  delete Rhs;
  delete mftfp;
  Vlevel1("**** DarcySNES destructor **** (end)",verbose);
}

static DarcySNES* static_ds_ptr = 0;

void
DarcySNES::SetTheDarcySNES(DarcySNES* ptr) 
{
    static_ds_ptr = ptr;
}

int 
DarcySNES::Solve(PArray<MultiFab>& RhoSat_new,
		 PArray<MultiFab>& RhoSat_old,
		 int               rhosatComp,
		 PArray<MultiFab>& P_new,
		 int               pComp,
		 Real              delta_t)
{
  Vlevel1("**** DarcySNES::Solve **** (begin)",verbose);

  // Set up data structures for state variables
  ResetRhoSat(RhoSat_new,RhoSat_old,rhosatComp,P_new,pComp);

  MFTower& RhsMFT = *Rhs;
  MFTower& SolnMFT = *Pnew;
  dt = delta_t;

  // Copy from MFTowers in state to Vec structures
  PetscErrorCode ierr;
  ierr = layout.MFTowerToVec(RhsV,RhsMFT,0); CHKPETSC(ierr);
  ierr = layout.MFTowerToVec(SolnV,SolnMFT,0); CHKPETSC(ierr);

  DarcySNES::SetTheDarcySNES(this);
  ierr = SNESSolve(snes,PETSC_NULL,SolnV); CHKPETSC(ierr);
  DarcySNES::SetTheDarcySNES(0);

  int iters;
  ierr = SNESGetIterationNumber(snes,&iters);CHKPETSC(ierr);

  SNESConvergedReason reason;
  ierr = SNESGetConvergedReason(snes,&reason);

  if (reason <= 0) {
    Vlevel1("**** DarcySNES::Solve **** (end)",verbose);
    if (abort_on_nl_fail) {
      BoxLib::Abort();
    }
    return reason;
  }

  // Copy solution from Vec back into state
  ierr = layout.VecToMFTower(SolnMFT,SolnV,0); CHKPETSC(ierr);
  
  Vlevel1("**** DarcySNES::Solve **** (end)",verbose);
  return reason;
}

void
DarcySNES::ResetRhoSat(PArray<MultiFab>& S_new,
                       PArray<MultiFab>& S_old,
		       int               rsComp,
                       PArray<MultiFab>& P_new,
		       int               pComp)
{
  Vlevel2("**** DarcySNES::ResetRhoSat **** (begin)",verbose);
  delete RhoSatOld; RhoSatOld = new MFTower(layout,S_old,nLevs,rsComp,1);
  delete RhoSatNew; RhoSatNew = new MFTower(layout,S_new,nLevs,rsComp,1);
  delete Pnew; Pnew = new MFTower(layout,P_new,nLevs,pComp,1);
  rhoSat_comp = 0;
  pressure_comp = 0; // This is offset within MFTower by its BaseComp
  Vlevel2("**** DarcySNES::ResetRhoSat **** (end)",verbose);
}

void
DarcySNES::BuildOpSkel(Mat& J)
{
  Vlevel3("**** DarcySNES::BuildOpSkel **** (begin)",verbose);
  int num_rows = 1;
  int rows[1]; // At the moment, only set one row at a time
  Array<Real> vals;
  Array<int> cols;
  
  const Array<BoxArray>& gridArray = layout.GridArray();
  const Array<IntVect>& refRatio = layout.RefRatio();
  const PArray<Layout::MultiNodeFab>& nodes = layout.Nodes();
  const PArray<Layout::MultiIntFab>& nodeIds = layout.NodeIds();
  const Array<Array<IVSMap> >& growCellStencil = mftfp->GrowCellStencil();
  int nLevs = layout.NumLevels();
  
  PetscErrorCode ierr;
  Layout::IntFab reg_neighbors;
  std::set<int> neighbors;
  typedef BaseFab<std::set<int> > ISetFab;
  typedef FabArray<ISetFab> MultiSetFab;
  PArray<MultiSetFab> crseContribs(nLevs,PArrayManage);
  
  int myproc = ParallelDescriptor::MyProc();
  int numprocs = ParallelDescriptor::NProcs();
  
  for (int lev=nLevs-1; lev>=0; --lev) 
    {
      const Array<IVSMap>& growCellStencilLev = growCellStencil[lev];
      const Layout::MultiNodeFab& nodeLev = nodes[lev];
      const Layout::MultiIntFab& nodeIdsLev = nodeIds[lev];

      Layout::MultiIntFab crseIds; // coarse cell ids at fine grid, distributed per fine patches
      crseContribs.set(lev,new MultiSetFab);
      if (lev>0) {
	BoxArray bacg = BoxArray(gridArray[lev]).coarsen(refRatio[lev-1]).grow(1);
	crseIds.define(bacg,1,0,Fab_allocate);
            
	const Layout::MultiIntFab& crseIds_orig = nodeIds[lev-1]; // crse cells through periodic boundary
	BoxArray gcba = BoxArray(crseIds_orig.boxArray()).grow(crseIds_orig.nGrow());
	Layout::MultiIntFab tmp(gcba,1,0);
	for (MFIter mfi(crseIds_orig); mfi.isValid(); ++mfi) {
	  tmp[mfi].copy(crseIds_orig[mfi]); // NOTE: Assumes grow cells already filled
	}
	crseIds.copy(tmp); // Parallel copy

	crseContribs[lev].define(bacg,1,0,Fab_allocate);
      }

      std::map<IntVect,std::set<int>,IntVect::Compare> stencil;
      if (lev<nLevs-1) {
	// Pack up the crseContribs for a parallel copy
	const BoxArray& ba = gridArray[lev];
	MultiSetFab& crseContribsFine = crseContribs[lev+1];
        const DistributionMapping& dm = nodeLev.DistributionMap();
	std::map<int,Array<int> > ccArrays;
	for (MFIter mfi(crseContribsFine); mfi.isValid(); ++mfi) {
	  const ISetFab& ccFab = crseContribsFine[mfi];
	  const Box& vbox = mfi.validbox();
	  std::vector< std::pair<int,Box> > isects = ba.intersections(vbox);
	  for (int i=0; i<isects.size(); ++i) {
            int dst_proc = dm[isects[i].first];

            // HACK  This was originally written for parallel, but when I tried it in serial, the entire 
            // crseContribs structure was ignored!!  For now, set this up as a communication, even if 
            // serial...probably an easy logic issue to clear up....famous last words...
	    if (1 || dst_proc != myproc) {
	      for (IntVect iv(vbox.smallEnd()), iEnd=vbox.bigEnd(); iv<=iEnd; vbox.next(iv))
		{
		  const std::set<int>& ids = ccFab(iv,0);
		  int thisSize = ids.size();
		  if (thisSize) {
		    Array<int>& ints = ccArrays[dst_proc];
		    int old_cc_size = ints.size();
		    int delta_cc = BL_SPACEDIM + 1 + ids.size();
		    int new_cc_size = old_cc_size + delta_cc;

		    ints.resize(new_cc_size);
		    for (int d=0; d<BL_SPACEDIM; ++d) {
		      ints[old_cc_size+d] = iv[d];
		    }
		    ints[old_cc_size+BL_SPACEDIM] = ids.size();
		    int cnt=0;
		    for (std::set<int>::const_iterator it=ids.begin(), End=ids.end(); it!=End; ++it, ++cnt) {
		      ints[old_cc_size+BL_SPACEDIM+1+cnt] = *it;
		    }
		  }
		}
	    }
	  }
	}

	int total_num_to_send = 0;
	Array<int> sends(numprocs,0);
	Array<int> soffsets(numprocs,0);
	for (int i=0; i<numprocs; ++i) {
	  sends[i] = ccArrays[i].size();
	  total_num_to_send += sends[i];
	  if (i>0) {
	    soffsets[i] = soffsets[i-1] + ccArrays[i-1].size();
	  }
	}
	Array<int> sbuf(total_num_to_send);
	for (int i=0; i<numprocs; ++i) {
	  for (int j=0; j<ccArrays[i].size(); ++j) {
	    sbuf[soffsets[i] + j] = ccArrays[i][j];
	  }
	}

	Array<int> recvs(numprocs);
	BL_MPI_REQUIRE( MPI_Alltoall(sends.dataPtr(),
				     1,
				     ParallelDescriptor::Mpi_typemap<int>::type(),
				     recvs.dataPtr(),
				     1,
				     ParallelDescriptor::Mpi_typemap<int>::type(),
				     ParallelDescriptor::Communicator()) );
            
	int total_num_to_recv = 0;
	Array<int> roffsets(numprocs,0);
	for (int i=0; i<numprocs; ++i) {
	  total_num_to_recv += recvs[i];
	  if (i>0) {
	    roffsets[i] = roffsets[i-1] + recvs[i-1];
	  }
	}
	Array<int> rbuf(total_num_to_recv);
	BL_MPI_REQUIRE( MPI_Alltoallv(total_num_to_send == 0 ? 0 : sbuf.dataPtr(),
				      sends.dataPtr(),
				      soffsets.dataPtr(),
				      ParallelDescriptor::Mpi_typemap<int>::type(),
				      total_num_to_recv == 0 ? 0 : rbuf.dataPtr(),
				      recvs.dataPtr(),
				      roffsets.dataPtr(),
				      ParallelDescriptor::Mpi_typemap<int>::type(),
				      ParallelDescriptor::Communicator()) );
            
	for (int i=0; i<numprocs; ++i) {
	  int jcnt = roffsets[i];
	  while (jcnt < roffsets[i] + recvs[i]) {
	    IntVect iv(&(rbuf[jcnt]));
	    int size = rbuf[jcnt+BL_SPACEDIM];
	    std::set<int>& iset = stencil[iv];
	    for (int k=0; k<size; ++k) {
	      iset.insert(rbuf[jcnt+BL_SPACEDIM+1+k]);
	    }
	    jcnt += BL_SPACEDIM+1+size;
	  }
	}
      }

      for (MFIter mfi(nodeLev); mfi.isValid(); ++mfi) {
	const Layout::NodeFab& nodeFab = nodeLev[mfi];
	const Layout::IntFab& nodeIdFab = nodeIdsLev[mfi];
	const Layout::IntFab* crseIdFab = (lev>0  ?  &(crseIds[mfi])  : 0);
	const Box& vbox = mfi.validbox();

	for (IntVect iv(vbox.smallEnd()), iEnd=vbox.bigEnd(); iv<=iEnd; vbox.next(iv))
	  {
	    const Node& nC = nodeFab(iv,0);
	    if (nC.type==Node::VALID) {
	      rows[0] = nodeIdFab(iv,0);
	      neighbors.clear();

	      std::map<IntVect,std::set<int>,IntVect::Compare>::const_iterator sit=stencil.find(iv);
	      if (sit!=stencil.end()) {
		const std::set<int>& iset = sit->second;
		neighbors.insert(iset.begin(),iset.end());
	      }
	      neighbors.insert(rows[0]);

	      for (int d=0; d<BL_SPACEDIM; ++d) {
		for (int pm = -1; pm<2; pm+=2) {
		  std::set<int> nd;
		  IntVect ivA = iv  +  pm * BoxLib::BASISV(d);
		  IVScit it=growCellStencilLev[d].find(ivA);
		  if (it!=growCellStencilLev[d].end()) {
		    const Stencil& s = it->second;
		    for (Stencil::const_iterator it=s.begin(), End=s.end(); it!=End; ++it) {
		      const Node& node = it->first;
		      const IntVect& ivs = node.iv;
		      int slev = node.level;
		      if (slev==lev) {
			BL_ASSERT(nodeIdFab.box().contains(ivs));
			int idx = nodeIdFab(ivs,0);
			if (ivs != iv && idx>=0) { // idx<0 is Dirichlet data, iv added above
			  nd.insert(idx);
			}
		      }
		      else if (slev==lev-1) {
			BL_ASSERT(crseIdFab);
			BL_ASSERT(crseIdFab->box().contains(ivs));
			nd.insert((*crseIdFab)(ivs,0));
		      }
		      else {
			std::cout << "stencil: " << s << std::endl;
			BoxLib::Abort("Bad stencil");
		      }
		    }

		    // contribute to coarse cell stencil, if appropriate
		    const Node& offcenter_node = nodeFab(ivA,0);
		    if (offcenter_node.type==Node::VALID  &&  offcenter_node.level==lev-1) {
		      crseContribs[lev][mfi](offcenter_node.iv,0).insert(rows[0]);
		      crseContribs[lev][mfi](offcenter_node.iv,0).insert(nd.begin(),nd.end());
		    }
		  }
		  else {
		    int idx = nodeIdFab(ivA,0);
		    if (idx>=0) { // idx<0 is a covered cell
		      neighbors.insert(idx);
		    }
		  }

		  // Merge this arm into full set
		  neighbors.insert(nd.begin(),nd.end());

		}
	      }

	      int num_cols = -1;
              num_cols = neighbors.size();
              cols.resize(num_cols);
              vals.resize(num_cols,0);
              int cnt = 0;
              for (std::set<int>::const_iterator it=neighbors.begin(), End=neighbors.end(); it!=End; ++it) {
                cols[cnt++] = *it;
              }

	      ierr = MatSetValues(J,num_rows,rows,num_cols,cols.dataPtr(),vals.dataPtr(),INSERT_VALUES); CHKPETSC(ierr);
	    }
	  }
      }
    }

  ierr = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
  ierr = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY); CHKPETSC(ierr);
  Vlevel3("**** DarcySNES::BuildOpSkel **** (end)",verbose);
}

void
DarcySNES::FillPatch(MFTower& mft,
                     int sComp,
                     int nComp,
                     bool do_piecewise_constant)
{
  Vlevel3("**** DarcySNES::FillPatch **** (begin)",verbose);
  mftfp->FillGrowCells(mft,sComp,nComp,do_piecewise_constant,nLevs);
  Vlevel3("**** DarcySNES::FillPatch **** (end)",verbose);
}

Real Kr_given_Seff_Mualem(Real se, Real m, Real mInv) 
{
  Real Kr;
  if (se<=0) {
    Kr = 0;
  }
  else if (se>=1) {
    Kr = 1;
  }
  else {
    Real tmp = (1 - std::pow(1 - std::pow(se,mInv), m));
    Kr = std::sqrt(se)*tmp*tmp;
    Kr = std::min(1., std::max(0., Kr) );
  }
  return Kr;
}

Real Se_given_P_vanGenuchten(Real pressure, Real n, Real m, Real sigma)
{
  Real pcap = - pressure;
  Real se = (pcap > 0  ? std::pow(1 + std::pow(pcap*sigma,n),-m) : 1);
  return se;
}

void
DarcySNES::ReducedSaturationGivenPressure(const FArrayBox&      pressure,
					  int                   pComp,
					  FArrayBox&            reduced_saturation,
					  int                   satComp,
					  const Box&            box,
					  const Layout::IntFab& mask,
                                          Real                  maskedVal)
{
  Vlevel4("**** DarcySNES::ReducedSaturationGivenPressure **** (begin)",verbose);
  FArrayBox tp(box,1); tp.copy(pressure,pComp,0,1);
  Layout::IntFab ti(box,1); ti.copy(mask);
  FArrayBox ts(box,1);
  const Real* pdat = tp.dataPtr();
  Real* satdat = ts.dataPtr();
  const int* mdat = ti.dataPtr();
  int nPts = box.numPts();
  Real n = 1/(1-m);
  for (int i=0; i<nPts; ++i) {
    satdat[i] = ( mdat[i]>=0  ?  Se_given_P_vanGenuchten(pdat[i],n,m,sigma)  : maskedVal);
  }
  reduced_saturation.copy(ts,0,satComp,1);
  Vlevel4("**** DarcySNES::ReducedSaturationGivenPressure **** (end)",verbose);
}

void
DarcySNES::KrGivenReducedSaturation(const FArrayBox& reduced_saturation,
                                    int              satComp,
                                    FArrayBox&       Kr,
                                    int              kComp,
                                    const Box&       box,
                                    const Layout::IntFab& mask,
                                    Real                  maskedVal)
{
  Vlevel4("**** DarcySNES::KrGivenReducedSaturation **** (begin)",verbose);
  FArrayBox ts(box,1); ts.copy(reduced_saturation,satComp,0,1);
  Layout::IntFab ti(box,1); ti.copy(mask);
  FArrayBox tk(box,1);
  const Real* sdat = ts.dataPtr();
  Real* kdat = tk.dataPtr();
  const int* mdat = ti.dataPtr();
  int nPts = box.numPts();

  Real mInv = 1/m;

  Real DKrDse_sThresh = 0;
  Real DKrDse_tot = 0;
  Real Dse_tot = 1 - se_threshold;
  if (se_threshold < 1 && se_threshold>0) {
    Real Kr_se_threshold = Kr_given_Seff_Mualem(se_threshold,m,mInv);
    Real se_eps = 1.e-3;
    Real sep = std::min(1.,se_threshold + se_eps);
    Real Kr_se_thresholdp = Kr_given_Seff_Mualem(sep,m,mInv);
    DKrDse_sThresh = -(Kr_se_thresholdp - Kr_se_threshold)/(sep-se_threshold);
    DKrDse_tot = (Kr_se_threshold - 1)/Dse_tot;
  }
  
  for (int i=0; i<nPts; ++i) {
    if (mdat[i]>=0) {
      if (sdat[i] >= 1) {
        kdat[i] = 1;
      }
      else if (sdat[i] <= se_threshold) { 
        kdat[i] = (sdat[i] <= 0 ?  0  :  Kr_given_Seff_Mualem(sdat[i],m,mInv));
      }
      else {
        Real Dse = 1 - sdat[i];
        kdat[i] = 1 + Dse*Dse * DKrDse_tot / Dse_tot
          + Dse*Dse * (Dse-Dse_tot) * (DKrDse_sThresh-2*DKrDse_tot)/(Dse_tot*Dse_tot);
      }
    }
    else {
      kdat[i] = maskedVal;
    }
  }
  Kr.copy(tk,0,kComp,1);
  Vlevel4("**** DarcySNES::KrGivenReducedSaturation **** (end)",verbose);
}

void
DarcySNES::RhoSatGivenReducedSaturation(const FArrayBox&      reduced_saturation,
                                        int                   redComp,
                                        FArrayBox&            RhoSat,
                                        int                   rsComp,
                                        const Box&            box,
                                        const Layout::IntFab& mask,
                                        Real                  maskedVal)
{
  Vlevel4("**** DarcySNES::RhoSatGivenReducedSaturation **** (begin)",verbose);
  FArrayBox tred(box,1); tred.copy(reduced_saturation,redComp,0,1);
  Layout::IntFab ti(box,1); ti.copy(mask);
  FArrayBox trs(box,1);

  const Real* reddat = tred.dataPtr();
  Real* rsdat = trs.dataPtr();
  const int* mdat = ti.dataPtr();
  int nPts = box.numPts();

  Real oneMinusSr = 1 - Sr;
  for (int i=0; i<nPts; ++i) {
    Real se = reddat[i];
    if (mdat[i]>=0) {
      if (se<=0) {
        rsdat[i] = density * Sr;
      }
      else if (se>=1) {
        rsdat[i] = density;
      }
      else {
        rsdat[i] = density * ( Sr +  se * oneMinusSr );
      }
    }
    else {
      rsdat[i] = maskedVal;
    }
  }
  RhoSat.copy(trs,0,rsComp,1);
  Vlevel4("**** DarcySNES::RhoSatGivenReducedSaturation **** (end)",verbose);
}

void
DarcySNES::ComputeDarcyFlux(PArray<MFTower>& darcy_flux,
                            int              fComp,
                            MFTower&         pressure,
                            int              pComp,
                            MFTower&         rhoSat,
                            int              rsComp)
{
  Vlevel3("**** DarcySNES::ComputeDarcyFlux **** (begin)",verbose);
  int nGrowOp = 1;
  BL_ASSERT(pressure.NGrow() >= nGrowOp);
  mlb.SetDirichletValues(pressure,pComp);

  // Convert grow cells of pressure into extrapolated values
  bool do_piecewise_constant = false;
  int nComp = 1;
  FillPatch(pressure,pComp,nComp,do_piecewise_constant);

  // Note: For computing se(p), rhosat(se,rho) and Kr(s), we make use of a "mask"
  //       to allow us to easily avoid computing with invalid data,
  //       ie, "covered" data that may not have workable values
  //       This mask is built easily with the help of the layouts 
  //       ability to provide node numbers.  For cells inside the
  //       domain at least, valid cells will have idx>=0.  Outside
  //       the domain, we assume that "SetDirichletValues" has 
  //       put/left something computable.

  Layout::IntFab ifab;
  for (int lev=0; lev<nLevs; ++lev) {
    const MultiFab& pmf = pressure[lev];
    MultiFab& smf = rhoSat[lev];
    for (MFIter mfi(pmf); mfi.isValid(); ++mfi) {
      const FArrayBox& pfab = pmf[mfi];
      FArrayBox& sfab = smf[mfi];
      const Box& gbox = BoxLib::grow(mfi.validbox(),nGrowOp);
      ifab.resize(gbox,1);
      ifab.setVal(0);
      const Box& obox = gbox & layout.GeomArray()[lev].Domain();
      layout.SetNodeIds(ifab,lev,mfi.index(),obox);
      const Box& ovlp = sfab.box() & gbox;
      ReducedSaturationGivenPressure(pfab,pComp+pressure.BaseComp(),sfab,rsComp+rhoSat.BaseComp(),ovlp,ifab); // For the moment rhoSat holds se
    }
  }

  // Make sure fine-fine cells have correct periodically-mapped values
  for (int lev=0; lev<nLevs; ++lev) {
    const Geometry& geom = layout.GeomArray()[lev];

    MultiFab& smf = rhoSat[lev];
    smf.FillBoundary(rsComp+rhoSat.BaseComp(),nComp, geom.periodicity());

    MultiFab& pmf = pressure[lev];
    pmf.FillBoundary(pComp+pressure.BaseComp(),nComp, geom.periodicity());
  }


  FArrayBox kr_tmp;
  Array<FArrayBox*> flux(BL_SPACEDIM);
  int fc[BL_SPACEDIM];
  for (int lev=0; lev<nLevs; ++lev) {
    const MultiFab& pmf = pressure[lev];
    MultiFab& smf = rhoSat[lev];

    const Geometry& geom = layout.GeomArray()[lev];
    const Real* dx = geom.CellSize();

    for (MFIter mfi(pmf); mfi.isValid(); ++mfi) {
      const FArrayBox& pfab = pmf[mfi];
      FArrayBox& sfab = smf[mfi];
      for (int d=0; d<BL_SPACEDIM; ++d) {
        flux[d] = &(darcy_flux[d][lev][mfi]);
        fc[d] = darcy_flux[d].BaseComp() + fComp;
      }
      const Box& vbox = mfi.validbox();
      const Box& gbox = BoxLib::grow(vbox,nGrowOp);
      ifab.resize(gbox,1);
      ifab.setVal(0);
      const Box& obox = gbox & geom.Domain();
      layout.SetNodeIds(ifab,lev,mfi.index(),obox);
      int iComp = 0;

      kr_tmp.resize(gbox,1);
      int krComp = 0;
      KrGivenReducedSaturation(sfab, rsComp+rhoSat.BaseComp(), kr_tmp, krComp, gbox, ifab);
      RhoSatGivenReducedSaturation(sfab, rsComp+rhoSat.BaseComp(), sfab, rsComp+rhoSat.BaseComp(), gbox, ifab);

      FORT_DARCYFLUX(pfab.dataPtr(pComp+pressure.BaseComp()),ARLIM(pfab.loVect()), ARLIM(pfab.hiVect()),
		     kr_tmp.dataPtr(krComp),ARLIM(kr_tmp.loVect()), ARLIM(kr_tmp.hiVect()),
		     ifab.dataPtr(iComp),ARLIM(ifab.loVect()), ARLIM(ifab.hiVect()),
		     D_DECL(flux[0]->dataPtr(fc[0]),flux[1]->dataPtr(fc[1]),flux[2]->dataPtr(fc[2])),
		     D_DECL(ARLIM(flux[0]->loVect()),ARLIM(flux[1]->loVect()),ARLIM(flux[2]->loVect())),
		     D_DECL(ARLIM(flux[0]->hiVect()),ARLIM(flux[1]->hiVect()),ARLIM(flux[2]->hiVect())),
		     vbox.loVect(), vbox.hiVect(), &kappax, &kappaz, &density, &mu, gravity.dataPtr(), dx);
    }
  }

  // Overwrite fluxes at boundary with boundary conditions
  mlb.SetInflowFlux(darcy_flux,fComp);

  // Average down fluxes
  if (nLevs>1) {
    int nComp = 1;
    for (int d=0; d<BL_SPACEDIM; ++d) {
      MFTower::AverageDown(darcy_flux[d],fComp,nComp,nLevs);
    }    
  }
  Vlevel3("**** DarcySNES::ComputeDarcyFlux **** (end)",verbose);
}

void
DarcySNES::DivRhoU(MFTower& DivRhoU,
		   int      druComp,
                   MFTower& pressure,
		   int      pComp,
                   MFTower& rhoSat,
                   int      rsComp)
{
  Vlevel3("**** DarcySNES::DivRhoU **** (begin)",verbose);
  int fComp = 0;

  ComputeDarcyFlux(DarcyFlux,fComp,pressure,pComp,rhoSat,rsComp);

  // Get the divergence of the Darcy Flux
  int mult = 1;
  int nComp = 1;
  MFTower::ECtoCCdiv(DivRhoU,DarcyFlux,mult,fComp,druComp,nComp,nLevs);

  Vlevel3("**** DarcySNES::DivRhoU **** (end)",verbose);
}


void
DarcySNES::DpDtResidual(MFTower& residual,
			int      resComp,
                        MFTower& pressure,
			int      pComp,
                        Real     dt)
{
  Vlevel2("**** DarcySNES::DpDtResidual **** (begin)",verbose);

  DivRhoU(residual,resComp,pressure,pComp,*RhoSatNew,rhoSat_comp);

  int nComp=1;

  for (int lev=0; lev<nLevs; ++lev)
  {
      MultiFab& Rlev = residual[lev];
      for (MFIter mfi(Rlev); mfi.isValid(); ++mfi) {
	const Box& vbox = mfi.validbox();
	FArrayBox& Res = Rlev[mfi];
	const FArrayBox& RSn = (*RhoSatOld)[lev][mfi];
	const FArrayBox& RSnp1 = (*RhoSatNew)[lev][mfi];
	FORT_RS_PDOTRES(Res.dataPtr(resComp+residual.BaseComp()),ARLIM(Res.loVect()), ARLIM(Res.hiVect()),
			RSn.dataPtr(RhoSatOld->BaseComp()),ARLIM(RSn.loVect()), ARLIM(RSn.hiVect()),
			RSnp1.dataPtr(RhoSatNew->BaseComp()),ARLIM(RSnp1.loVect()), ARLIM(RSnp1.hiVect()),
		        &dt, &porosity, vbox.loVect(), vbox.hiVect(), &nComp);
      }
  }
  Vlevel2("**** DarcySNES::DpDtResidual **** (end)",verbose);
}

Real TotalVolume()
{
  const RealBox& rb = Geometry::ProbDomain();
  Real vol = 1;
  for (int d=0; d<BL_SPACEDIM; ++d) {
    vol *= rb.length(d);
  }
  return vol;
}

PetscErrorCode 
DarcyRes_DpDt(SNES snes,Vec x,Vec f,void *dummy)
{
  PetscErrorCode ierr; 
  DarcySNES* ds = static_cast<DarcySNES*>(dummy);
  if (!ds) {
    BoxLib::Abort("Bad cast in DarcyRes_DpDt");
  }
  
  MFTower& xMFT = *(ds->Pnew);
  MFTower& fMFT = *(ds->Rhs);
  
  Layout& layout = ds->layout;
  Real dt=ds->dt;
  ierr = layout.VecToMFTower(xMFT,x,ds->pressure_comp); CHKPETSC(ierr);
  
  int resComp = 0;
  ds->DpDtResidual(fMFT,resComp,xMFT,ds->pressure_comp,dt);

#if 1
    // Scale residual by cell volume/sqrt(total volume)
  Real sqrt_total_volume_inv = std::sqrt(1/TotalVolume());
  int nComp = 1;
  int nLevs = ds->nLevs;
  for (int lev=0; lev<nLevs; ++lev)
  {
    MultiFab::Multiply(fMFT[lev],layout.Volume(lev),0,0,nComp,0);
    fMFT[lev].mult(sqrt_total_volume_inv,0,1);
  }
#endif
  
  ierr = layout.MFTowerToVec(f,fMFT,ds->pressure_comp); CHKPETSC(ierr);
  
  return 0;
}



#if defined(PETSC_3_2)
#include <private/matimpl.h>
#else
#include <petsc-private/matimpl.h> 
#endif

PetscErrorCode  
DarcyMatFDColoringApply(Mat J,MatFDColoring coloring,Vec x1,MatStructure *flag,void *sctx)
{
  Vlevel3("**** DarcyMatFDColoringApply **** (begin)",verbose_SNES);
  PetscErrorCode (*f)(void*,Vec,Vec,void*) = (PetscErrorCode (*)(void*,Vec,Vec,void *))coloring->f;
  PetscErrorCode ierr;
  PetscInt       k,start,end,l,row,col,srow,m1,m2;
  PetscScalar    *y,*w3_array,*w4_array,*w5_array;
  PetscReal      epsilon = coloring->error_rel;
  Vec            w1=coloring->w1,w2=coloring->w2,w3;
  void           *fctx = coloring->fctx;
  PetscBool      flg = PETSC_FALSE;
  PetscInt       ctype=coloring->ctype,N;
  Vec            x1_tmp;

  PetscFunctionBegin;    
  PetscValidHeaderSpecific(J,MAT_CLASSID,1);
  PetscValidHeaderSpecific(coloring,MAT_FDCOLORING_CLASSID,2);
  PetscValidHeaderSpecific(x1,VEC_CLASSID,3);
  if (!f) SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_WRONGSTATE,"Must call MatFDColoringSetFunction()");

  ierr = PetscLogEventBegin(MAT_FDColoringApply,coloring,J,x1,0);CHKPETSC(ierr);
  ierr = MatSetUnfactored(J);CHKPETSC(ierr);
  ierr = PetscOptionsGetBool(PETSC_NULL,"-mat_fd_coloring_dont_rezero",&flg,PETSC_NULL);CHKPETSC(ierr);
  if (flg) {
    ierr = PetscInfo(coloring,"Not calling MatZeroEntries()\n");CHKPETSC(ierr);
  } else {
    PetscBool  assembled;
    ierr = MatAssembled(J,&assembled);CHKPETSC(ierr);
    if (assembled) {
      ierr = MatZeroEntries(J);CHKPETSC(ierr);
    }
  }

  x1_tmp = x1; 
  if (!coloring->vscale){ 
    ierr = VecDuplicate(x1_tmp,&coloring->vscale);CHKPETSC(ierr);
  }
    
  /*
    This is a horrible, horrible, hack. See DMMGComputeJacobian_Multigrid() it inproperly sets
    coloring->F for the coarser grids from the finest
  */
  if (coloring->F) {
    ierr = VecGetLocalSize(coloring->F,&m1);CHKPETSC(ierr);
    ierr = VecGetLocalSize(w1,&m2);CHKPETSC(ierr);
    if (m1 != m2) {  
      coloring->F = 0; 
      }    
    }   


  DarcySNES* ds = static_ds_ptr;
  BL_ASSERT(ds);

  ierr = VecGetOwnershipRange(w1,&start,&end);CHKPETSC(ierr); /* OwnershipRange is used by ghosted x! */
      
  /* Set w1 = F(x1) */
  if (coloring->F) {
    w1          = coloring->F; /* use already computed value of function */
    coloring->F = 0; 
  } else {
    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
    ierr = (*f)(sctx,x1_tmp,w1,fctx);CHKPETSC(ierr);
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
  }
      
  if (!coloring->w3) {
    ierr = VecDuplicate(x1_tmp,&coloring->w3);CHKPETSC(ierr);
    ierr = PetscLogObjectParent(coloring,coloring->w3);CHKPETSC(ierr);
  }
  w3 = coloring->w3;

    /* Compute all the local scale factors, including ghost points */
  ierr = VecGetLocalSize(x1_tmp,&N);CHKPETSC(ierr);

  ierr = VecSet(coloring->vscale,1/epsilon);

  if (ctype == IS_COLORING_GLOBAL){
      ierr = VecGhostUpdateBegin(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKPETSC(ierr);
      ierr = VecGhostUpdateEnd(coloring->vscale,INSERT_VALUES,SCATTER_FORWARD);CHKPETSC(ierr);
  }
  
  if (coloring->vscaleforrow) {
  } else SETERRQ(((PetscObject)J)->comm,PETSC_ERR_ARG_NULL,"Null Object: coloring->vscaleforrow");

  /*
    Loop over each color, the perturbation is a simple constant, epsilon
  */
  Vec w4, w5;
  if (centered_diff_J) {
    w4 = ds->Wrk1V;
    w5 = ds->Wrk2V;
  }

  for (k=0; k<coloring->ncolors; k++) { 
    coloring->currentcolor = k;
    
    ierr = VecCopy(x1_tmp,w3);CHKPETSC(ierr);
    if (centered_diff_J) {
      ierr = VecCopy(x1_tmp,w4);CHKPETSC(ierr);
      ierr = VecCopy(x1_tmp,w5);CHKPETSC(ierr);
    }

    ierr = VecGetArray(w3,&w3_array);CHKPETSC(ierr);
    if (centered_diff_J) {
      ierr = VecGetArray(w4,&w4_array);CHKPETSC(ierr);
      ierr = VecGetArray(w5,&w5_array);CHKPETSC(ierr);
    }

    if (ctype == IS_COLORING_GLOBAL) {
      w3_array = w3_array - start;
      if (centered_diff_J) {
        w4_array = w4_array - start;
        w5_array = w5_array - start;
      }
    }

    for (l=0; l<coloring->ncolumns[k]; l++) {
      col = coloring->columns[k][l];    /* local column of the matrix we are probing for */
      w3_array[col] += epsilon;   // w3 = x1 + dx
      if (centered_diff_J) {
        w4_array[col] -= epsilon; // w4 = x1 - dx
      } 
    } 

    if (ctype == IS_COLORING_GLOBAL) {
      w3_array = w3_array + start;
      if (centered_diff_J) {
        w4_array = w4_array + start;
        w5_array = w5_array + start;
      }
    }
    ierr = VecRestoreArray(w3,&w3_array);CHKPETSC(ierr);
    if (centered_diff_J) {
      ierr = VecRestoreArray(w4,&w4_array);CHKPETSC(ierr);
      ierr = VecRestoreArray(w5,&w5_array);CHKPETSC(ierr);
    }

    ierr = PetscLogEventBegin(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);
    ierr = (*f)(sctx,w3,w2,fctx);CHKPETSC(ierr);                       // w2 = F(w3) = F(x1 + dx)
    if (centered_diff_J) {ierr = (*f)(sctx,w4,w5,fctx);CHKPETSC(ierr);}  // w5 = F(w4) = F(x1 - dx)
    ierr = PetscLogEventEnd(MAT_FDColoringFunction,0,0,0,0);CHKPETSC(ierr);

    PetscReal epsilon_inv = 1/epsilon;
    if (centered_diff_J) {
      epsilon_inv *= 0.5;
      ierr = VecAXPY(w2,-1.0,w5);CHKPETSC(ierr); // w2 = F(x1 + dx) - F(x1 - dx)
    } else {
      ierr = VecAXPY(w2,-1.0,w1);CHKPETSC(ierr); // w2 = F(x1 + dx) - F(x1)
    }
    
    // Insert (w2_j / dx) into J_ij
    ierr = VecGetArray(w2,&y);CHKPETSC(ierr);          
    for (l=0; l<coloring->nrows[k]; l++) {
      row    = coloring->rows[k][l];             /* local row index */
      col    = coloring->columnsforrow[k][l];    /* global column index */
      y[row] *= epsilon_inv;                     /* dx = epsilon */
      srow   = row + start;                      /* global row index */
      ierr   = MatSetValues(J,1,&srow,1,&col,y+row,INSERT_VALUES);CHKPETSC(ierr);
    }
    ierr = VecRestoreArray(w2,&y);CHKPETSC(ierr);
    
  } /* endof for each color */
  
   
  coloring->currentcolor = -1;
  ierr  = MatAssemblyBegin(J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  ierr  = MatAssemblyEnd(J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  ierr = PetscLogEventEnd(MAT_FDColoringApply,coloring,J,x1,0);CHKPETSC(ierr);

  Vlevel3("**** DarcyMatFDColoringApply **** (end)",verbose_SNES);
  PetscFunctionReturn(0);
}

PetscErrorCode 
DarcyComputeJacobianColor(SNES snes,Vec x1,Mat *J,Mat *B,MatStructure *flag,void *ctx)
{
  Vlevel3("**** DarcyComputeJacobianColor **** (begin)",verbose_SNES);

  MatFDColoring  color = (MatFDColoring) ctx;
  PetscErrorCode ierr;
  Vec            f;
  PetscErrorCode (*ff)(void),(*fd)(void);

  PetscFunctionBegin;
  PetscValidHeaderSpecific(color,MAT_FDCOLORING_CLASSID,6);

  *flag = SAME_NONZERO_PATTERN;
  ierr  = SNESGetFunction(snes,&f,(PetscErrorCode (**)(SNES,Vec,Vec,void*))&ff,0);CHKPETSC(ierr);
  ierr  = MatFDColoringGetFunction(color,&fd,PETSC_NULL);CHKPETSC(ierr);
  if (fd == ff) { /* reuse function value computed in SNES */
    ierr  = MatFDColoringSetF(color,f);CHKPETSC(ierr);
  }
  ierr = DarcyMatFDColoringApply(*B,color,x1,flag,snes);CHKPETSC(ierr);
  if (*J != *B) {
    ierr = MatAssemblyBegin(*J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
    ierr = MatAssemblyEnd(*J,MAT_FINAL_ASSEMBLY);CHKPETSC(ierr);
  }

  Vlevel3("**** DarcyComputeJacobianColor **** (end)",verbose_SNES);
  PetscFunctionReturn(0);
}
