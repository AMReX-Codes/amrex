#include <iostream.h>
#include <fstream.h>
#include <math.h>
#include <stdio.h>

#ifdef BL_ARCH_CRAY
#ifdef BL_USE_DOUBLE
#error "DOUBLE PRECISION NOT ALLOWED ON CRAY"
#endif
#endif

#ifndef        WIN32
#include <unistd.h>
#endif

#include <Utility.H>
#include <Tracer.H>
#include <Box.H>
#include <BoxArray.H>
#include <Geometry.H>
#include <ParmParse.H>

#include <TestMCViscBndry.H>
#include <DivVis.H>
#include <LO_BCTYPES.H>
#include <MCMultiGrid.H>
#include <MCCGSolver.H>

#ifdef BL_USE_NEW_HFILES
#include <new>
using std::setprecision;
#ifndef WIN32
using std::set_new_handler;
#endif
#else
#include <new.h>
#endif

const int NPROCS = 1;

void inspectFAB( FARRAYBOX& unfab )
{
    char filename[256];
    sprintf(filename,"inspect.fab");
    ofstream osout(filename,ios::out);
    unfab.writeOn(osout);
    osout.close();
#if BL_SPACEDIM == 2
    system("/usr/local/bin/amrvis2d -fab inspect.fab&");
#else
    system("/usr/local/bin/amrvis3d -fab inspect.fab&");
#endif    
}

void
ViewMultiFab(ostream &os, const MultiFab& phi,
	     const REAL* H, const BoxList& bs, const BOX& container,
	     int ratio=2, REAL bgValue=0.0);

void inspectMF( MultiFab& soln )
{
    ofstream os("pltfile");
    int comp = 0;
    REAL s_min=soln.min(comp);
    REAL s_max=soln.max(comp);
    REAL bgValue = s_min - 0.1*(s_max - s_min);
    FArrayBox::setFormat(FABio::FAB_NATIVE);
    REAL h[BL_SPACEDIM];
    h[0] = h[1] = 1.;
#if BL_SPACEDIM>2
    h[2] = 1.;
#endif
    const BoxArray &ba = soln.boxArray();
    Box container = ba.minimalBox();
    int ratio = 1;
    const BoxList bl(ba);
    ViewMultiFab(os,soln,h,bl,container,ratio,bgValue);
}

#include <main_F.H>

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
REAL* fabdat = (fab).dataPtr();

BoxList
readBoxList(aString file, BOX& domain );

main(int   argc,
     char* argv[])
{
  //
  // Make sure to catch new failures.
  //
#ifndef WIN32
  set_new_handler(Utility::OutOfMemory);
#endif
  
  TRACER("mg");
  
  if(argc < 2) {
      cerr << "usage:  " << argv[0] << " inputsfile [options]" << '\n';
      exit(-1);
  }
  
  // Parse command line
  ParmParse pp(argc-2,argv+2,NULL,argv[1]); 
  
  int nprocs = NPROCS; pp.query("nprocs", nprocs);
#ifndef BL_USE_BSP
  if (nprocs > 1)
  {
      cerr << "Error in main:  multiple processors specified with "
           << "code compiled without a parallel library.\n";
      exit(-1);
  }
#endif
  ParallelDescriptor::StartParallel(nprocs);
  
  int n;

#if BL_SPACEDIM == 2
  Box container(IntVect(0,0),IntVect(11,11));
  aString boxfile("gr.2_small_a") ;
#elif BL_SPACEDIM == 3
  Box container(IntVect(0,0,0),IntVect(11,11,11));
  aString boxfile("grids/gr.3_small_a") ;
#endif
  pp.query("boxes", boxfile);
  BoxArray bs(readBoxList(boxfile,container));
  Geometry geom(container);
  REAL H[BL_SPACEDIM];
  for (n=0; n<BL_SPACEDIM; n++) {
    H[n] = ( geom.ProbHi(n) - geom.ProbLo(n) )/(container.length(n));
  } // -->> over dimension
  int ratio=2; pp.query("ratio", ratio);

  // allocate/init soln and rhs
  int Ncomp=BL_SPACEDIM;
  int Nghost=0;
  int Ngrids=bs.length();
  MultiFab soln(bs, Ncomp, Nghost, Fab_allocate); soln.setVal(0.0);
  MultiFab out(bs, Ncomp, Nghost, Fab_allocate); 
  MultiFab rhs(bs, Ncomp, Nghost, Fab_allocate); 
  rhs.setVal(0.0);
  for(MultiFabIterator rhsmfi(rhs); rhsmfi.isValid(); ++rhsmfi) {
    DEF_LIMITS(rhsmfi(),fdat,flo,fhi);
    FORT_FILLRHS(fdat,ARLIM(flo),ARLIM(fhi),H,&Ncomp);
  }

  // Create the boundary object
  MCViscBndry vbd(bs,geom);

  // Create the BCRec's interpreted by ViscBndry objects
#if BL_SPACEDIM==2
  Array<BCRec> pbcarray(4);
  pbcarray[0] = BCRec(D_DECL(REFLECT_ODD,REFLECT_EVEN,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[1] = BCRec(D_DECL(REFLECT_EVEN,REFLECT_ODD,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[2] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[3] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#elif BL_SPACEDIM==3
  Array<BCRec> pbcarray(12);
#define TEST_EXT 0
#if TEST_EXT
  pbcarray[0] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
  pbcarray[1] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
  pbcarray[2] = BCRec(EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR,EXT_DIR);
#else
  pbcarray[0] = BCRec(INT_DIR,INT_DIR,INT_DIR, INT_DIR,INT_DIR,INT_DIR);
  pbcarray[1] = BCRec(INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR);
  pbcarray[2] = BCRec(INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR,INT_DIR);
#endif
  // the other bc's don't really matter.  Just init with anything
  pbcarray[3] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[4] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[5] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[6] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[7] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[8] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[9] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[10] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
  pbcarray[11] = BCRec(D_DECL(EXT_DIR,EXT_DIR,EXT_DIR),
			  D_DECL(EXT_DIR,EXT_DIR,EXT_DIR));
#endif

  Nghost = 1; // need space for bc info
  MultiFab fine(bs,Ncomp,Nghost,Fab_allocate);
  for(MultiFabIterator finemfi(fine); finemfi.isValid(); ++finemfi) {
    DEF_LIMITS(finemfi(),fdat,flo,fhi);
    FORT_FILLFINE(fdat,ARLIM(flo),ARLIM(fhi),H,&Ncomp);
  }

  // Create "background coarse data"
  BOX crse_bx(container);
  crse_bx.coarsen(ratio);
  crse_bx.grow(1);
  REAL h_crse[BL_SPACEDIM];
  for (n=0; n<BL_SPACEDIM; n++) h_crse[n] = H[n]*ratio;
  FARRAYBOX crse_fab(crse_bx,Ncomp);
  DEF_LIMITS(crse_fab,cdat,clo,chi);
  FORT_FILLCRSE(cdat,ARLIM(clo),ARLIM(chi),h_crse,&Ncomp);
  
  // Create coarse boundary register, fill w/data from coarse FAB
  int bndry_InRad=0;
  int bndry_OutRad=1;
  int bndry_Extent=1;
  BoxArray cbs(bs); cbs.coarsen(ratio);
  BndryRegister cbr(cbs,bndry_InRad,bndry_OutRad,bndry_Extent,Ncomp);
  for (OrientationIter face; face; ++face) {
    Orientation f = face();
    FabSet& bnd_fs(cbr[f]);
    bnd_fs.copyFrom(crse_fab);
  }
  
  // Interpolate crse data to fine boundary, where applicable
  int cbr_Nstart=0;
  int fine_Nstart=0;
  int bndry_Nstart=0;
  vbd.setBndryValues(cbr,cbr_Nstart,fine,fine_Nstart,
		     bndry_Nstart,Ncomp,ratio,pbcarray);
  
  Nghost = 1; // other variables don't need extra space

  DivVis lp(bs,vbd,H);

  REAL a = 0.0;
  REAL b[BL_SPACEDIM];
  b[0] = 1.0;
  b[1] = 1.0;
#if BL_SPACEDIM>2
  b[2] = 1.0;
#endif
  MultiFab  acoefs;
  acoefs.define(bs, Ncomp, Nghost, Fab_allocate);
  acoefs.setVal(a);
  MultiFab bcoefs[BL_SPACEDIM];
  for (n=0; n<BL_SPACEDIM; ++n) {
    BoxArray bsC(bs);
    bcoefs[n].define(bsC.surroundingNodes(n), 1,
		     Nghost, Fab_allocate);
#if 1
    for(MultiFabIterator bmfi(bcoefs[n]); bmfi.isValid(); ++bmfi) {
      DEF_LIMITS(bmfi(),fdat,flo,fhi);
      FORT_MAKEMU(fdat,ARLIM(flo),ARLIM(fhi),H,n);
    }
#else
    bcoefs[n].setVal(b[n]);
#endif
  } // -->> over dimension
  lp.setCoefficients(acoefs, bcoefs);
#if 1
  lp.maxOrder(4);
#endif

  Nghost = 1;
  MultiFab tsoln(bs, Ncomp, Nghost, Fab_allocate); 
  tsoln.setVal(0.0);
#if 1
  tsoln.copy(fine);
#endif
#if 1
  // testing apply
  lp.apply(out,tsoln);
  Box subbox = out[0].box();
  REAL n1 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
  cout << "n1 output is "<<n1<<endl;
  out.minus(rhs,0,BL_SPACEDIM,0);
  // special to single grid prob
  REAL n2 = out[0].norm(subbox,1,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
  cout << "n2 difference is "<<n2<<endl;
#if 0
  subbox.grow(-1);
  REAL n3 = out[0].norm(subbox,0,0,BL_SPACEDIM)*pow(H[0],BL_SPACEDIM);
  cout << "n3 difference is "<<n3<<endl;
#endif
  ofstream outos("dogfile");
  FArrayBox::setFormat(FABio::FAB_NATIVE);
  REAL tmin = out.min(0);
  tmin -= 1.0;
  const BoxList bl(bs);
  ViewMultiFab(outos,out, H, bl, container, ratio, tmin);
  outos.close();
#endif
#if 1
  // testing flux computation
  BoxArray xfluxbox(bs);
  xfluxbox.surroundingNodes(0);
  MultiFab xflux(xfluxbox,Ncomp,Nghost,Fab_allocate);
  xflux.setVal(1.e30);
  BoxArray yfluxbox(bs);
  yfluxbox.surroundingNodes(1);
  MultiFab yflux(yfluxbox,Ncomp,Nghost,Fab_allocate);
  yflux.setVal(1.e30);
#if BL_SPACEDIM>2
  BoxArray zfluxbox(bs);
  zfluxbox.surroundingNodes(2);
  MultiFab zflux(zfluxbox,Ncomp,Nghost,Fab_allocate);
  zflux.setVal(1.e30);
#endif
  lp.compFlux(xflux,
	      yflux,
#if BL_SPACEDIM>2
	      zflux,
#endif
	      tsoln);

  MultiFab &flux = yflux;
  ofstream fluxos("fluxfile");
  REAL fmin = flux.min(0);
  fmin -= 1.0;
  Box fcontainer(container);
  fcontainer.grow(1);
  const BoxList flxbl(flux.boxArray());
  ViewMultiFab(fluxos,flux, H, flxbl, fcontainer, ratio, fmin);
  fluxos.close();
#endif



  ofstream ofs;

  REAL tolerance = 1.0e-10; pp.query("tol", tolerance);
  REAL tolerance_abs = 1.0e-10; pp.query("tol_abs", tolerance_abs);

#if 0
  bool use_mg_pre = false;
  MCCGSolver cg(lp,use_mg_pre);
  cg.solve(soln,rhs,tolerance,tolerance_abs);
#else
  MCMultiGrid mg(lp);
  mg.solve(soln,rhs,tolerance,tolerance_abs);
#endif

#if 0
  // apply operator to soln to see if really satisfies eqn
  tsoln.copy(soln);
  lp.apply(out,tsoln);
  soln.copy(out);
#endif
# if 0
  // output a big fab instead of a plt file
  Box shrunkcon = container;
  shrunkcon.grow(-1);
  FArrayBox motherfab(shrunkcon,BL_SPACEDIM);
  motherfab.setVal(0.);
  soln.copy(motherfab);
  ofstream os("dog.fab");
  motherfab.writeOn(os);
  os.close();
#elif 0
  // Output multifab directly
  ofstream os("dog.mfab");
  soln.writeOn(os);
  os.close();
#else
  for(MultiFabIterator smfi(soln); smfi.isValid(); ++smfi) {
      DependentMultiFabIterator fmfi(smfi,fine);
      smfi() -= fmfi();
  }
  for( int icomp=0; icomp < BL_SPACEDIM ; icomp++ ){
    cout << icomp << "  "<<soln.min(icomp) << " " << soln.max(icomp)<<endl;
  }

  ofstream os("pltfile");
  REAL bg_val=1.0;
  FArrayBox::setFormat(FABio::FAB_NATIVE);
  ViewMultiFab(os,soln, H, bl, container, ratio, bg_val);
  os.close();
#endif

  EndParallel();
    
}

BoxList
readBoxList(const aString file, BOX& domain )
{
    BoxList retval;
    ifstream boxspec(file.c_str());
    if( !boxspec ) {
	BoxLib::Error("readBoxList: unable to open " + *file.c_str());
    }
    boxspec >> domain;
    
    int numbox;
    boxspec >> numbox;

    for(int i=0; i<numbox; i++) {
	BOX tmpbox;
	boxspec >> tmpbox;
	if( ! domain.contains(tmpbox)) {
	    cerr << "readBoxList: bogus box " << tmpbox << endl;
	    exit(1);
        }
	retval.append(tmpbox);
    }
    boxspec.close();
    return retval;
}


void
ViewMultiFab(ostream &os, const MultiFab& phi,
	     const REAL* H, const BoxList& bs, const BOX& container,
	     int ratio, REAL bgValue)
{
    int n;
    int nvar = phi.nComp();
    Tuple< REAL, BL_SPACEDIM > probLo;
    Tuple< REAL, BL_SPACEDIM > probHi;
    const int* lo=container.loVect();
    const int* hi=container.hiVect();
    REAL HBG[BL_SPACEDIM];
    for (n=0; n<BL_SPACEDIM; n++) {
	probLo[n]=H[n]*lo[n];
	probHi[n]=H[n]*(hi[n]+1.0);
	HBG[n]=H[n]*ratio;
    }
    BOX bxbg(container);
    bxbg.coarsen(ratio);
    FARRAYBOX fab(bxbg,nvar);
    fab.setVal(bgValue);


    os << nvar         << "\n"; // Number of states dumped
    for(int ivar=0; ivar<nvar; ivar++){
      char dog[256];
      sprintf(dog,"dummy%d",ivar);
      os << dog << "\n";
    }
    os << BL_SPACEDIM     << "\n"; // Dimension of data
    os << "0.0"        << "\n"; // Simulation time of dump
    os << "1"          << "\n"; // Finest level dumped
    
    for (n=0; n<BL_SPACEDIM; n++)  os << probLo[n]  << " ";
    os << "\n";                                  // Position of lo-sides
    for (n=0; n<BL_SPACEDIM; n++)  os << probHi[n]  << " "; 
    os << "\n";                                  // Position of hi-sides
    os << ratio << "\n";                         // 0:f_lev-2 ratio
    os << bxbg << " " << container << "\n";      // 0:f_lev-1 container
    os << "0 0\n";                               // 0:f_lev-1 nSteps
    for (n=0; n<BL_SPACEDIM; n++) os << HBG[n] << " ";
    os << "\n";                                  // Grid spacing on background
    for (n=0; n<BL_SPACEDIM; n++) os << H[n]   << " "; 
    os << "\n";                                  // Grid spacing of data
    os << "0\n";                                 // Coord sys flag (0=cart)
    os << "0\n";                                 // BC flag (0=no BC info)

      // dump base grid
    os << "0 1\n";                      // level and number of grids
    os << bxbg << "\n";                 // For each grid, dump box
    os << "0\n";                        //                level
    os << "0\n";                        //                steps
    os << "0.0\n";                      //                time
    for (n=0; n<BL_SPACEDIM; n++) { 
        os << probLo[n] << " "          //                probLo
           << probHi[n] << "\n";        //                probHi
    }                                   
    fab.writeOn(os,0,nvar);             //                fab dump valid data

      // dump actual data
    int ngrds=bs.length();
    os << "1 " << ngrds << "\n";        // level and number of grids
    const BoxArray bxa(bs);
    int i;
    for (i = 0; i<ngrds; i++) {
        const BOX& grd=bxa[i];          //
        os << grd << "\n";              // For each grid, dump box
	os << "1\n";                    //                level
	os << "0\n";                    //                steps
	os << "0.0\n";                  //                time
	for (n=0; n<BL_SPACEDIM; n++) {
	    os << probLo[n] << " "      //                probLo
               << probHi[n] << "\n";    //                probHi
	} // endfor: dimension

	FARRAYBOX dat(grd,nvar);
	dat.copy(phi[i],0,0,nvar);
	dat.writeOn(os,0,nvar);           //                 fab dump valid data
    }

}

