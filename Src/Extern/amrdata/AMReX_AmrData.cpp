// ---------------------------------------------------------------
// AmrData.cpp
// ---------------------------------------------------------------

#include <AMReX_AmrData.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_VisMF.H>
#include <AMReX_MFCopyDescriptor.H>

#include <iostream>
#include <string>
using std::string;
using std::cout;
using std::cerr;
using std::endl;
using std::strlen;

#include <iostream>
#include <fstream>
#include <cstdio>
using std::ios;
using std::ifstream;

#ifdef SHOWVAL
#undef SHOWVAL
#endif

#define SHOWVAL(val) { cout << #val << " = " << val << endl; }

#ifdef VSHOWVAL
#undef VSHOWVAL
#endif

#define VSHOWVAL(verbose, val) { if(verbose) { \
                 cout << #val << " = " << val << endl; } }


#if defined( BL_FORT_USE_UPPERCASE )
#  if (BL_SPACEDIM == 1)
#    define   FORT_PCINTERP    PCINTERP1D
#  elif (BL_SPACEDIM == 2)
#    define   FORT_CINTERP     CINTERP2D
#    define   FORT_PCINTERP    PCINTERP2D
#    define   FORT_CARTGRIDMINMAX CARTGRIDMINMAX2D
#  elif (BL_SPACEDIM == 3)
#    define   FORT_CINTERP     CINTERP3D
#    define   FORT_PCINTERP    PCINTERP3D
#    define   FORT_CARTGRIDMINMAX CARTGRIDMINMAX3D
#  endif
#elif defined( BL_FORT_USE_LOWERCASE )
#  if (BL_SPACEDIM == 1)
#    define   FORT_PCINTERP    pcinterp1d
#  elif (BL_SPACEDIM == 2)
#    define   FORT_CINTERP     cinterp2d
#    define   FORT_PCINTERP    pcinterp2d
#    define   FORT_CARTGRIDMINMAX cartgridminmax2d
#  elif (BL_SPACEDIM == 3)
#    define   FORT_CINTERP     cinterp3d
#    define   FORT_PCINTERP    pcinterp3d
#    define   FORT_CARTGRIDMINMAX cartgridminmax3d
#  endif
#else
#  if (BL_SPACEDIM == 1)
#    define   FORT_PCINTERP    pcinterp1d_
#  elif (BL_SPACEDIM == 2)
#    define   FORT_CINTERP     cinterp2d_
#    define   FORT_PCINTERP    pcinterp2d_
#    define   FORT_CARTGRIDMINMAX cartgridminmax2d_
#  elif (BL_SPACEDIM == 3)
#    define   FORT_CINTERP     cinterp3d_
#    define   FORT_PCINTERP    pcinterp3d_
#    define   FORT_CARTGRIDMINMAX cartgridminmax3d_
#  endif
#endif


extern "C" {
#if (BL_SPACEDIM != 1)
  void FORT_CINTERP(amrex::Real *fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                  const int *fblo, const int *fbhi,
                  const int &nvar, const int &lratio,
		  const amrex::Real *crse, const int &clo, const int &chi,
		  const int *cslo, const int *cshi,
		  const int *fslo, const int *fshi,
		  amrex::Real *cslope, const int &c_len,
		  amrex::Real *fslope, amrex::Real *fdat, const int &f_len,
		  amrex::Real *foff);
#endif

  void FORT_PCINTERP(amrex::Real *fine, AMREX_ARLIM_P(flo), AMREX_ARLIM_P(fhi),
                   const int *fblo, const int *fbhi,
		   const int &lrat, const int &nvar,
		   const amrex::Real *crse, AMREX_ARLIM_P(clo), AMREX_ARLIM_P(chi),
		   const int *cblo, const int *cbhi,
		   amrex::Real *temp, const int &tlo, const int &thi);

#if (BL_SPACEDIM != 1)
  void FORT_CARTGRIDMINMAX (amrex::Real *data, AMREX_ARLIM_P(dlo), AMREX_ARLIM_P(dhi),
		            const amrex::Real *vfrac, const amrex::Real &vfeps,
		            amrex::Real &dmin, amrex::Real &dmax);
#endif
}

namespace amrex {

bool AmrData::verbose = false;
int  AmrData::skipPltLines  = 0;
int  AmrData::sBoundaryWidth = 0;

// ---------------------------------------------------------------
AmrData::AmrData() {
  probSize.resize(BL_SPACEDIM, -1.0);
  probLo.resize(BL_SPACEDIM,  0.0);
  probHi.resize(BL_SPACEDIM, -1.0);
  plotVars.clear();
  nRegions = 0;
  boundaryWidth = 0;
}


// ---------------------------------------------------------------
AmrData::~AmrData() {
   for(int lev(0); lev < regions.size(); ++lev) {
     for(int i(0); i < regions[lev].size(); ++i) {
       delete regions[lev][i];
     }
   }

   for(int lev(0); lev < dataGrids.size(); ++lev) {
     for(int i(0); i< dataGrids[lev].size(); ++i) {
       delete dataGrids[lev][i];
     }
   }

   for(int lev(0); lev < visMF.size(); ++lev) {
     for(int i(0); i < visMF[lev].size(); ++i) {
       delete visMF[lev][i];
     }
   }
}


// ---------------------------------------------------------------
namespace {
  void mytrim(char *str) {
    int i(strlen(str));
    for(int n(i - 1); n >= 0; --n) {
      if( str[n] > ' ' ) {
        break;
      }
      str[n] = 0;
    }
  }
}


// ---------------------------------------------------------------
bool AmrData::ReadData(const string &filename, Amrvis::FileType filetype) {
   fileType = filetype;
   bCartGrid = false;
   bShowBody = false;
   vCartGrid = -1;
   bTerrain = false;
   if(filetype == Amrvis::FAB || filetype == Amrvis::MULTIFAB) {
     return ReadNonPlotfileData(filename, filetype);
   }
   if(filetype == Amrvis::PROFDATA) {
      return false;    // ---- profdata will be handled later
   }

   int i, j, k, width;
   fileName = filename;

   MFInfo Fab_noallocate;
   Fab_noallocate.SetAlloc(false);

    string File = filename;

#ifdef BL_NO_PARALLEL_IO
    // do nothing
#else
    File += '/';
    File += "Header";
#endif

    VisMF::IO_Buffer io_buffer(VisMF::GetIOBufferSize());

    ifstream isPltIn;

    isPltIn.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

   if(verbose) {
     if(ParallelDescriptor::IOProcessor()) {
       cout << "AmrData::opening file = " << File << endl;
     }
   }

   isPltIn.open(File.c_str(), ios::in);
   if(isPltIn.fail()) {
     if(ParallelDescriptor::IOProcessor()) {
      cerr << "AmrData::Unable to open file:  " << File << endl;
     }
      return false;
   }

   char skipBuff[Amrvis::LINELENGTH];
   for(i = 0; i < skipPltLines; ++i) {
     isPltIn.getline(skipBuff, Amrvis::LINELENGTH);
     if(ParallelDescriptor::IOProcessor()) {
       cout << "Skipped line in pltfile = " << skipBuff << endl;
     }
   }

     isPltIn >> plotFileVersion;
     if(strncmp(plotFileVersion.c_str(), "CartGrid", 8) == 0) {
       bCartGrid = true;
       bShowBody = true;
     }
     if(strcmp(plotFileVersion.c_str(), "CartGrid-V1.0") == 0) {
       vCartGrid = 10;
     }
     if(strcmp(plotFileVersion.c_str(), "CartGrid-V1.2") == 0) {
       vCartGrid = 12;
     }
     if(strcmp(plotFileVersion.c_str(), "CartGrid-V2.0") == 0) {
       vCartGrid = 20;
     }
     if(bCartGrid && vCartGrid < 0) {
       if(ParallelDescriptor::IOProcessor()) {
         cerr << "**** Unknown CartGrid version:  " << plotFileVersion << endl;
       }
     }
     if(strncmp(plotFileVersion.c_str(), "Terrain", 7) == 0) {
       bTerrain = true;
     }
     if(verbose) {
       if(ParallelDescriptor::IOProcessor()) {
         cout << "Plot file version:  " << plotFileVersion << endl;
	 if(bCartGrid) {
	   cout << ":::: Found a CartGrid file:  version " << vCartGrid << endl;
	 }
	 if(bTerrain) {
	   cout << ":::: Found a Terrain file type." << endl;
	 }
       }
     }

     // read list of variables
     isPltIn >> nComp;
      if(nComp < 1 || nComp > 1024) {  // arbitrarily limit to 1024
        if(ParallelDescriptor::IOProcessor()) {
          cerr << "Error in AmrData:  bad nComp = " << nComp << endl;
        }
        return false;
      }
      if(ParallelDescriptor::IOProcessor()) {
        VSHOWVAL(verbose, nComp);
      }

      plotVars.resize(nComp);
      char plotVarName[Amrvis::LINELENGTH];
      bool bVFracFound(false);
      isPltIn.getline(plotVarName, Amrvis::LINELENGTH); // eat white space left by op<<
      for(i = 0; i < nComp; ++i) {
        isPltIn.getline(plotVarName, Amrvis::LINELENGTH);
	mytrim(plotVarName);
	if(bCartGrid) {  // check various permutations of vfrac name
	  if(strcmp(plotVarName, "vol_frac") == 0) {
	    cout << "+++++++++++++ found bad vfrac name:  " << plotVarName << endl;
	    strcpy(plotVarName, "vfrac");
	    cout << "+++++++++++++                  now:  " << plotVarName << endl;
	  }
	  if(strcmp(plotVarName, "volfrac") == 0) {
	    cout << "+++++++++++++ found bad vfrac name:  " << plotVarName << endl;
	    strcpy(plotVarName, "vfrac");
	    cout << "+++++++++++++                  now:  " << plotVarName << endl;
	  }
	  if(strcmp(plotVarName, "vfrac") == 0) {
            bVFracFound = true;
	  }
	}
        plotVars[i] = plotVarName;
        if(ParallelDescriptor::IOProcessor()) {
          VSHOWVAL(verbose, plotVarName);
	}
      }
      if(bCartGrid && bVFracFound == false) {
	cerr << endl << " ~~~~ Error:  no vfrac found for a " << plotFileVersion
	     << " file."  << endl << endl;
	return false;
      }

      int spacedim;
      isPltIn >>  spacedim >> time >> finestLevel;
      if(ParallelDescriptor::IOProcessor()) {
        VSHOWVAL(verbose, spacedim);
        VSHOWVAL(verbose, time);
        VSHOWVAL(verbose, finestLevel);
      }
      if(spacedim != BL_SPACEDIM) {
        if(ParallelDescriptor::IOProcessor()) {
	  cerr << endl << " ~~~~ Error:  You are using " << BL_SPACEDIM
	       << "D amrvis "
	       << "to look at a " << spacedim << "D file." << endl << endl;
	}
	return false;
      }
      if(finestLevel < 0) {
        if(ParallelDescriptor::IOProcessor()) {
          cerr << "Error in AmrData:  bad finestLevel = " << finestLevel << endl;
	}
        return false;
      }
      for(i = 0; i < BL_SPACEDIM; ++i) {
        isPltIn >> probLo[i];
        if(verbose) {
          if(ParallelDescriptor::IOProcessor()) {
	    cout << "probLo[" << i << "] = " << probLo[i] << endl;
	  }
	}
      }
      for(i = 0; i < BL_SPACEDIM; ++i) {
        isPltIn >> probHi[i];
        if(verbose) {
          if(ParallelDescriptor::IOProcessor()) {
	    cout << "probHi[" << i << "] = " << probHi[i] << endl;
	  }
	}
      }
      if(verbose) {
        if(ParallelDescriptor::IOProcessor()) {
	  if(finestLevel > 0) {
	    cout << "Resizing refRatio to size = " << finestLevel << endl;
	  }
	}
      }
      if(finestLevel == 0) {
        refRatio.resize(1, 1);
      } else {
        refRatio.resize(finestLevel, -1);
      }
      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }
      bool bIVRefRatio(false);
      if(isPltIn.peek() == '(') {  // it is an IntVect
        bIVRefRatio = true;
      }
      for(i = 0; i < finestLevel; ++i) {
	// try to guess if refRatio is an IntVect
	if(bIVRefRatio) {  // it is an IntVect
	  IntVect ivRefRatio;
	  isPltIn >> ivRefRatio;
          if(verbose) {
            if(ParallelDescriptor::IOProcessor()) {
	      cout << "IntVect refRatio[" << i << "] = " << ivRefRatio << endl;
	    }
	  }
	  refRatio[i] = ivRefRatio[0];  // non-uniform ref ratios not supported
	} else {
          isPltIn >> refRatio[i];
	}
        if(verbose) {
          if(ParallelDescriptor::IOProcessor()) {
	    cout << "refRatio[" << i << "] = " << refRatio[i] << endl;
	  }
	}
      }
      for(i = 0; i < finestLevel; ++i ) {
        if(refRatio[i] < 2 || refRatio[i] > 32 ) {
          if(ParallelDescriptor::IOProcessor()) {
            cerr << "Error in AmrData:  bad refRatio at level " << i << " = "
	         << refRatio[i] << endl;
	  }
          return false;
        }
      }
      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }
      probDomain.resize(finestLevel + 1);
      maxDomain.resize(finestLevel + 1);
      for(i = 0; i <= finestLevel; ++i) {
        isPltIn >> probDomain[i];
	if(verbose) {
          if(ParallelDescriptor::IOProcessor()) {
	    cout << "probDomain[" << i << "] = " << probDomain[i] << endl;
	  }
	}
        if( ! probDomain[i].ok()) {
          if(ParallelDescriptor::IOProcessor()) {
            cerr << "Error in AmrData:  bad probDomain[" << i << "] = "
	         << probDomain[i] << endl;
	  }
          return false;
        }
      }

      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }
#if 0      
      char lstepbuff[128];
      isPltIn.getline(lstepbuff, 128);  // ignore levelsteps--some files have
				   // finestlevel of these, others have
				   // finestlevel + 1

          
      if(verbose) {
        if(ParallelDescriptor::IOProcessor()) {
	  cout << "Ignored levelSteps = " << lstepbuff << endl;
	}
      }
#else
      // Get level steps - why do some have fl+1 ? (note above...just throw away last one here)
      levelSteps.resize(finestLevel+1);
      for(i = 0; i <= finestLevel; ++i) {
          isPltIn >> levelSteps[i];
      }
      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }

      if(verbose) {
        if(ParallelDescriptor::IOProcessor()) {
            for(i = 0; i <= finestLevel; ++i) {
                cout << "LevelSteps[" << i << "] = " << levelSteps[i] <<  endl;
            }
	}
      }
#endif
      
      dxLevel.resize(finestLevel + 1);
      for(i = 0; i <= finestLevel; ++i) {
        dxLevel[i].resize(BL_SPACEDIM);
        for(k = 0; k < BL_SPACEDIM; ++k) {
	  isPltIn >> dxLevel[i][k];
	  if(verbose) {
            if(ParallelDescriptor::IOProcessor()) {
	      cout << "dxLevel[" << i << "][" << k << "] = "
		   << dxLevel[i][k] << endl;
	    }
	  }
	}
      }

      vfEps.resize(finestLevel + 1);  // must resize these even if not cartGrid
      afEps.resize(finestLevel + 1);
      if(bCartGrid && vCartGrid < 20) {
        for(i = 0; i <= finestLevel; ++i) {
          isPltIn >> vfEps[i];
          if(verbose && ParallelDescriptor::IOProcessor()) {
            cout << "vfEps[" << i << "] = " << vfEps[i] << endl;
          }
        }
      }

      for(i = 0; i < BL_SPACEDIM; ++i) {
        probSize[i] = probHi[i] - probLo[i];
        if(probSize[i] <= 0.0 ) {
          if(ParallelDescriptor::IOProcessor()) {
            cerr << "Error in AmrData:  bad probSize[" << i << "] = "
	         << probSize[i] << endl;
	  }
          return false;
	}
      }

      isPltIn >> coordSys;
      if(ParallelDescriptor::IOProcessor()) {
        VSHOWVAL(verbose, coordSys);
      }
      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }

      isPltIn >> width;   // width of bndry regions
      if(ParallelDescriptor::IOProcessor()) {
        VSHOWVAL(verbose, width);
      }
      while(isPltIn.get() != '\n') {
        ;  // do nothing
      }

   dataGrids.resize(finestLevel + 1);
   dataGridsDefined.resize(finestLevel + 1);

   int lev;
   boundaryWidth = std::max(width, sBoundaryWidth);
   bool bRestrictDomain(maxDomain[0].ok());
   if(bRestrictDomain) {
      for(lev = 1; lev <= finestLevel; ++lev) {
        maxDomain[lev] = amrex::refine(maxDomain[lev-1],refRatio[lev-1]);
      }
   }
   Vector<Box> restrictDomain(finestLevel + 1);
   Vector<Box> extendRestrictDomain(finestLevel + 1);
   regions.resize(finestLevel + 1);
   for(lev = 0; lev <= finestLevel; ++lev) {
      restrictDomain[lev] = probDomain[lev];
      if(bRestrictDomain) {
        restrictDomain[lev] = maxDomain[lev];
      }
      extendRestrictDomain[lev] = amrex::grow(restrictDomain[lev],boundaryWidth);
      BoxList bndry_boxes = amrex::boxDiff(extendRestrictDomain[lev],
                                            restrictDomain[lev]);
      nRegions = bndry_boxes.size();

      BoxList::iterator bli = bndry_boxes.begin();
      regions[lev].resize(nRegions);
      i = 0;
      while(bli != bndry_boxes.end()) {
	regions[lev][i] = new FArrayBox(*bli, nComp);
	if(verbose) {
          if(ParallelDescriptor::IOProcessor()) {
	    cout << "BNDRY REGION " << i << " : " << *bli << endl;
	    cout << "    numPts = " << bli->numPts() << endl;
	  }
	}
	++i;
	++bli;
      }
   }

   // if positive set up and read bndry databoxes
   if(width > 0) {
     if(ParallelDescriptor::IOProcessor()) {
        cerr << "Error in AmrData:  Boundary width > 0 not supported:  width = "
	     << width << endl;
     }
     return false;
   }  // end if(width...)

   // read all grids but only save those inside the restricted region

    visMF.resize(finestLevel + 1);
    compIndexToVisMFMap.resize(nComp);
    compIndexToVisMFComponentMap.resize(nComp);
    gridLocLo.resize(finestLevel + 1);
    gridLocHi.resize(finestLevel + 1);

    for(i = 0; i <= finestLevel; ++i) {
      int nGrids;
      Real gTime;
      int iLevelSteps;
      isPltIn >> lev >> nGrids >> gTime >> iLevelSteps;
      if(ParallelDescriptor::IOProcessor()) {
        VSHOWVAL(verbose, lev);
        VSHOWVAL(verbose, nGrids);
        VSHOWVAL(verbose, gTime);
        VSHOWVAL(verbose, iLevelSteps);
      }
      if(i != lev) {
        if(ParallelDescriptor::IOProcessor()) {
	  cerr << "Level misrestart:mismatch on restart" << endl;
          cerr << "Error in AmrData:  Level mismatch:  read level " << lev
	       << " while expecting level " << i << endl;
	}
        return false;
      }
      if(nGrids < 1) {
        if(ParallelDescriptor::IOProcessor()) {
          cerr << "Error in AmrData:  bad nGrids = " << nGrids << endl;
	}
        return false;
      }

      gridLocLo[i].resize(nGrids);
      gridLocHi[i].resize(nGrids);
      for(int iloc = 0; iloc < nGrids; ++iloc) {
        gridLocLo[i][iloc].resize(BL_SPACEDIM);
        gridLocHi[i][iloc].resize(BL_SPACEDIM);
	for(int iDim = 0; iDim < BL_SPACEDIM; ++iDim) {
	  isPltIn >> gridLocLo[i][iloc][iDim] >>  gridLocHi[i][iloc][iDim];
          if(ParallelDescriptor::IOProcessor()) {
            VSHOWVAL(verbose, gridLocLo[i][iloc][iDim]);
            VSHOWVAL(verbose, gridLocHi[i][iloc][iDim]);
	  }
	}
      }

      // here we account for multiple multifabs in a plot file
      int currentIndexComp(0);
      int currentVisMF(0);
      dataGrids[i].resize(nComp);
      dataGridsDefined[i].resize(nComp);

      std::unique_ptr<DistributionMapping> dmap;

      while(currentIndexComp < nComp) {

        string mfNameRelative;
        isPltIn >> mfNameRelative;
        string mfName(fileName);
#ifdef BL_NO_PARALLEL_IO
        // do nothing
#else
        mfName += '/';
        mfName += mfNameRelative;
        VSHOWVAL(verbose, mfName);
        VSHOWVAL(verbose, mfNameRelative);
#endif

        visMF[i].resize(currentVisMF + 1);  // this preserves previous ones
        visMF[i][currentVisMF] = new VisMF(mfName);
	const BoxArray& ba = visMF[i][currentVisMF]->boxArray();
	if (!dmap) dmap.reset(new DistributionMapping(ba));
	int iComp(currentIndexComp);
        nGrow = visMF[i][currentVisMF]->nGrow();
        currentIndexComp += visMF[i][currentVisMF]->nComp();
	int currentVisMFComponent(0);
        for( ; iComp < currentIndexComp; ++iComp) {
          // make single component multifabs
          // defer reading the MultiFab data
	  dataGrids[i][iComp] = new MultiFab(ba, *dmap,
					     1,
			                     visMF[i][currentVisMF]->nGrow(),
					     Fab_noallocate);
          dataGridsDefined[i][iComp].resize(visMF[i][currentVisMF]->size(),
					    false);
          compIndexToVisMFMap[iComp] = currentVisMF;
          compIndexToVisMFComponentMap[iComp] = currentVisMFComponent;
          ++currentVisMFComponent;
        }

        ++currentVisMF;
      }  // end while

    }  // end for(i...finestLevel)

   // fill a set of temporary bndry regions surrounding the
   // restricted domain by extension from interior data
   // only use this data in bndry regions that did not
   // get better data from interior or input bndry regions
   for(lev = 0; lev <= finestLevel; ++lev) {
      Box inbox(restrictDomain[lev]);
      Box reg1(amrex::grow(restrictDomain[lev],boundaryWidth));
      Box reg2(amrex::grow(probDomain[lev],width));
      BoxList outside = amrex::boxDiff(reg1, reg2);
      if(outside.size() > 0) {
         // parts of the bndry have not been filled from the input
	 // data, must extending from interior regions

         for(int idir(0); idir < BL_SPACEDIM; ++idir) {
            Box bx(amrex::adjCellLo(inbox,idir,boundaryWidth));
	    Box br(bx);
	    for(k = 0; k < BL_SPACEDIM; ++k) {
	      if(k != idir) {
	        br.grow(k,boundaryWidth);
	      }
	    }

	    br.shift(idir,1);
	    FArrayBox tmpreg(br,nComp);
	    Box reg_bx = tmpreg.box();
	    reg_bx &= inbox;
	    FillInterior(tmpreg,lev,reg_bx);
	    br.shift(idir,-1);
	    FArrayBox tmpreg_lo(br,nComp);
	    tmpreg_lo.copy<RunOn::Host>(tmpreg,tmpreg.box(),0,tmpreg_lo.box(),0,nComp);

            // now fill out tmp region along idir direction
	    Box b_lo(amrex::adjCellLo(inbox,idir,1));
	    for(k = 1; k < boundaryWidth; ++k) {
	       Box btmp(b_lo);
	       btmp.shift(idir, -k);
	       tmpreg_lo.copy<RunOn::Host>(tmpreg_lo,b_lo,0,btmp,0,nComp);
	    }

	    // now fill out temp bndry region
	    Box b_src, b_dest;
	    int n;
	    for(k = 1; k < BL_SPACEDIM; ++k) {
	       int kdir = (idir + k) % BL_SPACEDIM;
 	       b_dest = amrex::adjCellLo(bx, kdir, 1);
	       b_src  = b_dest;
	       b_src  = b_src.shift(kdir, 1);
               for(n = 1; n <= boundaryWidth; ++n) {
	          tmpreg_lo.copy<RunOn::Host>(tmpreg_lo, b_src, 0, b_dest, 0, nComp);
	          b_dest.shift(kdir, -1);
	       }

	       b_dest = amrex::adjCellHi(bx,kdir,1);
	       b_src = b_dest;
	       b_src.shift(kdir,-1);
               for(n = 1; n <= boundaryWidth; ++n) {
	          tmpreg_lo.copy<RunOn::Host>(tmpreg_lo,b_src,0,b_dest,0,nComp);
	          b_dest.shift(kdir,1);
	       }
	       bx.grow(kdir,boundaryWidth);
	    }

	    // now copy into real bndry regions
	    for(j = 0; j < nRegions; ++j) {
	       FArrayBox *p = regions[lev][j];
	       Box p_box = p->box();
	       BoxList::iterator bli = outside.begin();
	       while(bli != outside.end()) {
                 Box ovlp(p_box);
		 ovlp &= *bli;
		 ovlp &= br;
		 if(ovlp.ok()) {
  		   p->copy<RunOn::Host>(tmpreg_lo, ovlp);
		 }
		 ++bli;
               }
	    }  // end for j

            // now work on the high side of the bndry region
            bx = amrex::adjCellHi(inbox,idir,boundaryWidth);
	    br = bx;
	    for(k = 0; k < BL_SPACEDIM; ++k) {
	      if(k != idir) br.grow(k, boundaryWidth);
	    }

	    br.shift(idir,-1);
	    FArrayBox tmpreg2(br,nComp);
	    reg_bx = tmpreg2.box();
	    reg_bx &= inbox;
	    FillInterior(tmpreg2,lev,reg_bx);
	    br.shift(idir,1);
	    FArrayBox tmpreg_hi(br,nComp);
	    tmpreg_hi.copy<RunOn::Host>(tmpreg2,tmpreg2.box(),0,tmpreg_hi.box(),0,nComp);

            // now fill out tmp region along idir direction
	    Box b_hi(amrex::adjCellHi(inbox,idir,1));
	    for(k = 1; k < boundaryWidth; ++k) {
	       Box btmp(b_hi);
	       btmp.shift(idir,k);
	       tmpreg_hi.copy<RunOn::Host>(tmpreg_hi,b_hi,0,btmp,0,nComp);
	    }

	    // now fill out temp bndry region
	    for(k = 1; k < BL_SPACEDIM; ++k) {
	       int kdir = (idir + k) % BL_SPACEDIM;
	       b_dest = amrex::adjCellLo(bx, kdir, 1);
	       b_src  = b_dest;
	       b_src.shift(kdir, 1);
               for(n = 1; n <= boundaryWidth; ++n) {
	          tmpreg_hi.copy<RunOn::Host>(tmpreg_hi, b_src, 0, b_dest, 0, nComp);
	          b_dest.shift(kdir,-1);
	       }

	       b_dest = amrex::adjCellHi(bx, kdir, 1);
	       b_src  = b_dest;
	       b_src.shift(kdir, -1);
               for(n = 1; n <= boundaryWidth; ++n) {
	          tmpreg_hi.copy<RunOn::Host>(tmpreg_hi, b_src, 0, b_dest, 0, nComp);
	          b_dest.shift(kdir, 1);
	       }
	       bx.grow(kdir, boundaryWidth);
	    }

	    // now copy into real bndry regions
	    for(j = 0; j < nRegions; ++j) {
	       FArrayBox *p = regions[lev][j];
	       Box p_box = p->box();
	       BoxList::iterator bli = outside.begin();
	       while(bli != outside.end()) {
                 Box ovlp(p_box);
		 ovlp &= *bli;
		 ovlp &= br;
		 if(ovlp.ok()) {
  		   p->copy<RunOn::Host>(tmpreg_hi, ovlp);
		 }
		 ++bli;
               }
	    }  // end for j

         }  // end for(idir...)
      }  // end if(outside.size())...

      outside.clear();

   }  // end for(lev...)

   if(bRestrictDomain) {
      Vector<Real> p_lo(BL_SPACEDIM), p_hi(BL_SPACEDIM);
      LoNodeLoc(0,maxDomain[0].smallEnd(),p_lo);
      HiNodeLoc(0,maxDomain[0].bigEnd(),p_hi);
      for(i = 0; i < BL_SPACEDIM; ++i) {
         probLo[i] = p_lo[i];
	 probHi[i] = p_hi[i];
	 probSize[i] = p_hi[i] - p_lo[i];
      }
      for(lev = 0; lev <= finestLevel; ++lev) {
         probDomain[lev] = maxDomain[lev];
      }
   }

   // ---- for CartGrid versions >= 20 write vfEps at the end of the header
   if(bCartGrid && vCartGrid >= 20) {
     for(i = 0; i <= finestLevel; ++i) {
       isPltIn >> vfEps[i];
       if(verbose && ParallelDescriptor::IOProcessor()) {
           cout << "vfEps[" << i << "] = " << vfEps[i] << endl;
       }
     }
   }
   isPltIn.close();

   return true;

}  // end ReadData


// ---------------------------------------------------------------
bool AmrData::ReadNonPlotfileData(const string &filename, Amrvis::FileType /*filetype*/) {
  const int LevelZero(0), LevelOne(1), BoxZero(0), ComponentZero(0);
  const int NVarZero(0), FabZero(0), IndexZero(0);
  const int iopNum(ParallelDescriptor::IOProcessorNumber());
  int i;
  if(verbose) {
    if(ParallelDescriptor::IOProcessor()) {
      cout << "AmrPlot::opening file = " << filename << endl;
    }
  }

  fileName = filename;

  time = 0;
  if(fileType == Amrvis::FAB) {
    finestLevel = LevelZero;
    plotFileVersion = "FromFAB";
  } else if(fileType == Amrvis::MULTIFAB) {
    finestLevel = 1;  // level zero is filler
    plotFileVersion = "FromMultiFAB";
  }
  probDomain.resize(finestLevel + 1);
  maxDomain.resize(finestLevel + 1);
  dxLevel.resize(finestLevel + 1);
  refRatio.resize(finestLevel + 1);
  if(fileType == Amrvis::FAB) {
    refRatio[LevelZero] = 1;
  } else if(fileType == Amrvis::MULTIFAB) {
    refRatio[LevelZero] = 2;
  }
  for(int iLevel(LevelZero); iLevel <= finestLevel; ++iLevel) {
    dxLevel[iLevel].resize(BL_SPACEDIM);
    for(i = 0; i < BL_SPACEDIM; ++i) {
      probLo[i] = 0.0;
      probHi[i] = 1.0;  // arbitrarily
      probSize[i] = probHi[i] - probLo[i];
      dxLevel[iLevel][i] = 0.0;  // temporarily
    }
  }

  dataGrids.resize(finestLevel + 1);
  dataGridsDefined.resize(finestLevel + 1);

  if(fileType == Amrvis::FAB) {
    ifstream is;
    is.open(filename.c_str(), ios::in);
    if(is.fail()) {
       cerr << "Unable to open plotfile: " << filename << endl;
       return false;
    }

    FArrayBox *newfab = new FArrayBox;
    nComp = newfab->readFrom(is, ComponentZero);  // read the first component
    Box fabbox(newfab->box());
    fabBoxArray.resize(1);
    fabBoxArray.set(BoxZero, fabbox);
    dataGrids[LevelZero].resize(nComp);
    dataGridsDefined[LevelZero].resize(nComp);
    dataGridsDefined[LevelZero][ComponentZero].resize(1);
    dataGrids[LevelZero][ComponentZero] = new MultiFab;
    int fabNGrow(0);

    Vector<int> pMap(fabBoxArray.size());
    pMap[BoxZero] = iopNum;
    DistributionMapping dMap(std::move(pMap));

    MFInfo Fab_noallocate;
    Fab_noallocate.SetAlloc(false);

    dataGrids[LevelZero][ComponentZero]->define(fabBoxArray, dMap, NVarZero,
                                                fabNGrow, Fab_noallocate);
    if(ParallelDescriptor::IOProcessor()) {
      dataGrids[LevelZero][ComponentZero]->setFab(FabZero, newfab);
    }
    dataGridsDefined[LevelZero][ComponentZero][IndexZero] = true;
    // read subsequent components
    // need to optimize this for lazy i/o
    for(int iComp = 1; iComp < nComp; ++iComp) {
      dataGrids[LevelZero][iComp] = new MultiFab;
      dataGrids[LevelZero][iComp]->define(fabBoxArray, dMap, NVarZero,
                                          fabNGrow, Fab_noallocate);
      if(ParallelDescriptor::IOProcessor()) {
        newfab = new FArrayBox;
        is.seekg(0, ios::beg);
        newfab->readFrom(is, iComp);  // read the iComp component
        dataGrids[LevelZero][iComp]->setFab(0, newfab);
      }
      dataGridsDefined[LevelZero][iComp].resize(1);
      dataGridsDefined[LevelZero][iComp][IndexZero] = true;
    }
    const int N(64);
    char fabname[N];  // arbitrarily
    plotVars.resize(nComp);
    for(i = 0; i < nComp; ++i) {
      if(snprintf(fabname, N, "%s%d", "Fab_", i) >= N) {
        amrex::Abort("AmrData::ReadNonPlotfileData: fabname buffer too small");
      }
      plotVars[i] = fabname;
    }
    probDomain[LevelZero] = newfab->box();
    for(i = 0; i < BL_SPACEDIM; ++i) {
      dxLevel[LevelZero][i] = 1.0 / probDomain[LevelZero].length(i);
    }
    is.close();

  } else if(fileType == Amrvis::MULTIFAB) {
    VisMF tempVisMF(filename);
    nComp = tempVisMF.nComp();
    probDomain[LevelOne] = tempVisMF.boxArray().minimalBox();
    probDomain[LevelZero] = probDomain[1];
    probDomain[LevelZero].coarsen(refRatio[0]);
    BoxArray mfBoxArray(tempVisMF.boxArray());
    BoxArray levelZeroBoxArray;
    levelZeroBoxArray.resize(1);
    levelZeroBoxArray.set(0, probDomain[LevelZero]);
    dataGrids[0].resize(nComp, NULL);
    dataGrids[1].resize(nComp, NULL);
    dataGridsDefined[LevelZero].resize(nComp);
    dataGridsDefined[LevelOne].resize(nComp);
    fabBoxArray.resize(1);
    fabBoxArray.set(BoxZero, probDomain[LevelZero]);

    int mfNGrow(0);
    const int N(64);
    char fabname[N];  // arbitrarily
    plotVars.resize(nComp);

    Vector<int> pMap(fabBoxArray.size());
    pMap[BoxZero] = iopNum;
    DistributionMapping dMap(std::move(pMap));

    MFInfo Fab_noallocate;
    Fab_noallocate.SetAlloc(false);

    for(int iComp(0); iComp < nComp; ++iComp) {
      if(snprintf(fabname, N, "%s%d", "MultiFab_", iComp) >= N) {
        amrex::Abort("AmrData::ReadNonPlotfileData: fabname buffer too small");
      }
      plotVars[iComp] = fabname;

      for(int iDim(0); iDim < BL_SPACEDIM; ++iDim) {
        dxLevel[0][iDim] = 1.0 / probDomain[0].length(iDim);
        dxLevel[1][iDim] = 1.0 / probDomain[1].length(iDim);
      }

      // set the level zero multifab
      dataGridsDefined[LevelZero][iComp].resize(1, false);
      dataGrids[LevelZero][iComp] = new MultiFab;
      dataGrids[LevelZero][iComp]->define(levelZeroBoxArray, dMap, NVarZero,
                                          mfNGrow, Fab_noallocate);
      if(ParallelDescriptor::IOProcessor()) {
        FArrayBox *newfab = new FArrayBox(probDomain[LevelZero], 1);
        Real levelZeroValue, zvMin, zvMax;
        zvMin = tempVisMF.min(0, iComp);  // init with first value
        zvMax = tempVisMF.max(0, iComp);  // init with first value
        for(int ic(0); ic < tempVisMF.size(); ++ic) {
	    zvMin = std::min(zvMin, tempVisMF.min(ic, iComp));
	    zvMax = std::max(zvMax, tempVisMF.max(ic, iComp));
        }
        levelZeroValue = zvMin;
        newfab->setVal<RunOn::Host>(levelZeroValue);
#if(BL_SPACEDIM == 2)
#ifdef BL_SETMFBACKGROUND
        Real *dptr = newfab->dataPtr();
        int idx;
        for(int icr(0); icr < newfab->box().length(1); ++icr) {
          for(int icc(0); icc < newfab->box().length(0); ++icc) {
	    idx = icc + (icr * newfab->box().length(0));
	    BL_ASSERT(idx < newfab->box().numPts());
	    if((icc + icr) % 5 == 0) {
              dptr[idx] = zvMax;
	    }
          }
        }
#endif
#endif
        dataGrids[LevelZero][iComp]->setFab(FabZero, newfab);
      }
      dataGridsDefined[LevelZero][iComp][IndexZero] = true;

    }  // end for(iComp...)


      // set the level one multifab

      // here we account for multiple multifabs in a plot file
      int currentIndexComp(0);
      int currentVisMF(0);
      visMF.resize(finestLevel + 1);
      compIndexToVisMFMap.resize(nComp);
      compIndexToVisMFComponentMap.resize(nComp);

      std::unique_ptr<DistributionMapping> dmap;

      while(currentIndexComp < nComp) {
        visMF[LevelOne].resize(currentVisMF + 1);  // this preserves previous ones
        visMF[LevelOne][currentVisMF] = new VisMF(filename);
	const BoxArray& ba = visMF[LevelOne][currentVisMF]->boxArray();
	if (!dmap) dmap.reset(new DistributionMapping(ba));
        int iComp(currentIndexComp);
        mfNGrow = visMF[LevelOne][currentVisMF]->nGrow();
        currentIndexComp += visMF[LevelOne][currentVisMF]->nComp();
        for(int currentVisMFComponent(0); iComp < currentIndexComp; ++iComp) {
          // make single component multifabs for level one
          dataGrids[1][iComp] =
	      new MultiFab(ba, *dmap,
			   1, visMF[LevelOne][currentVisMF]->nGrow(),
                           Fab_noallocate);
          dataGridsDefined[LevelOne][iComp].resize(visMF[LevelOne][currentVisMF]->size(), false);
          compIndexToVisMFMap[iComp] = currentVisMF;
          compIndexToVisMFComponentMap[iComp] = currentVisMFComponent;
          ++currentVisMFComponent;
        }

        ++currentVisMF;
      }  // end while

  }  // end if(fileType...)

  return true;
}


// ---------------------------------------------------------------
void AmrData::CellLoc(int lev, IntVect ix, Vector<Real> &pos) const {
   BL_ASSERT(pos.size() == dxLevel[lev].size());
   for(int i(0); i < BL_SPACEDIM; ++i) {
      pos[i] = probLo[i] + (dxLevel[lev][i])*(0.5 + Real(ix[i]));
   }
}


// ---------------------------------------------------------------
void AmrData::LoNodeLoc(int lev, IntVect ix, Vector<Real> &pos) const {
   BL_ASSERT(pos.size() == dxLevel[lev].size());
   for(int i(0); i < BL_SPACEDIM; ++i) {
      pos[i] = probLo[i] + (dxLevel[lev][i])*Real(ix[i]);
   }
}


// ---------------------------------------------------------------
void AmrData::HiNodeLoc(int lev, IntVect ix, Vector<Real> &pos) const {
   BL_ASSERT(pos.size() == dxLevel[lev].size());
   for(int i(0); i < BL_SPACEDIM; ++i) {
      pos[i] = probLo[i] + (dxLevel[lev][i])*Real(ix[i]+1);
   }
}


// ---------------------------------------------------------------
void AmrData::IntVectFromLocation(const int finestFillLevel,
                                  const Vector<Real> &location,
                                  IntVect &ivLoc, int &ivLevel,
				  IntVect &ivFinestFillLev)
{
   BL_ASSERT(location.size() == BL_SPACEDIM);
   BL_ASSERT(finestFillLevel <= finestLevel);

   int ival;

   for(int i(0); i < BL_SPACEDIM; ++i) {
      ival = probDomain[finestFillLevel].smallEnd()[i] +
                 (static_cast<int> ( (location[i] - probLo[i]) /
		                     dxLevel[finestFillLevel][i] ) );
      ivFinestFillLev.setVal(i, ival);
   }
   Box fflBox(ivFinestFillLev, ivFinestFillLev);
   ivLevel = FinestContainingLevel(fflBox, finestFillLevel);
   for(int i(0); i < BL_SPACEDIM; ++i) {
      ival = probDomain[ivLevel].smallEnd()[i] +
                 (static_cast<int> ( (location[i] - probLo[i]) /
		                     dxLevel[ivLevel][i] ) );
      ivLoc.setVal(i, ival);
   }
}


// ---------------------------------------------------------------
void AmrData::FillVar(FArrayBox *destFab, const Box &destBox,
		      int finestFillLevel, const string &varname, int procWithFabs)
{
  Vector<FArrayBox *> destFabs(1);
  Vector<Box> destBoxes(1);
  destFabs[0] = destFab;
  destBoxes[0] = destBox;

  FillVar(destFabs, destBoxes, finestFillLevel, varname, procWithFabs);
}


// ---------------------------------------------------------------
void AmrData::FillVar(MultiFab &destMultiFab, int finestFillLevel,
		      const string &varname, int destcomp)
{
  int numFillComps(1);
  Vector<string> varNames(numFillComps);
  Vector<int> destComps(numFillComps);
  varNames[0]  = varname;
  destComps[0] = destcomp;
  FillVar(destMultiFab, finestFillLevel, varNames, destComps);
}


// ---------------------------------------------------------------
void AmrData::FillVar(MultiFab &destMultiFab, int finestFillLevel,
		      const Vector<string> &varNames,
		      const Vector<int> &destFillComps)
{
// This function fills the destMultiFab which is defined on
// the finestFillLevel.

   BL_ASSERT(finestFillLevel >= 0 && finestFillLevel <= finestLevel);
   BoxArray destBoxes(destMultiFab.boxArray());
   for(int iIndex(0); iIndex < destBoxes.size(); ++iIndex) {
      if( ! probDomain[finestFillLevel].contains(destBoxes[iIndex])) {
         if(ParallelDescriptor::IOProcessor())  {
            cerr << "Error in AmrData::FillVar  -- probDomain does not contain destBoxes" << std::endl;
            cerr << "   Domain is  " << probDomain[finestFillLevel] << std::endl;
            cerr << "   ith box is " << destBoxes[iIndex]           << std::endl;
         }
         amrex::Abort("Error:  AmrData::FillVar");
      }
   }

    int myProc(ParallelDescriptor::MyProc());
    int srcComp(0);     // always 0 since AmrData uses single component MultiFabs
    int nFillComps(1);  // always
    int currentLevel;

    Vector<int> cumulativeRefRatios(finestFillLevel + 1, -1);

    cumulativeRefRatios[finestFillLevel] = 1;
    for(currentLevel = finestFillLevel - 1; currentLevel >= 0; --currentLevel) {
      cumulativeRefRatios[currentLevel] = cumulativeRefRatios[currentLevel + 1] *
                                          refRatio[currentLevel];
    }

    BL_ASSERT(varNames.size() == destFillComps.size());
    int nFillVars(varNames.size());

  for(int currentFillIndex(0); currentFillIndex < nFillVars; ++currentFillIndex) {
    int destComp(destFillComps[currentFillIndex]);
    int stateIndex(StateNumber(varNames[currentFillIndex]));
    // ensure the required grids are in memory
    for(currentLevel = 0; currentLevel <= finestFillLevel; ++currentLevel) {
      for(int iBox = 0; iBox < destBoxes.size(); ++iBox) {
	Box tempCoarseBox(destBoxes[iBox]);
        if(currentLevel != finestFillLevel) {
          tempCoarseBox.coarsen(cumulativeRefRatios[currentLevel]);
        }
        GetGrids(currentLevel, stateIndex, tempCoarseBox);
      }
    }

    MultiFabCopyDescriptor multiFabCopyDesc;
    Vector<MultiFabId> stateDataMFId(finestFillLevel + 1);
    for(currentLevel = 0; currentLevel <= finestFillLevel; ++currentLevel) {
      stateDataMFId[currentLevel] =
           multiFabCopyDesc.RegisterFabArray(dataGrids[currentLevel][stateIndex]);
    }

    BoxArray localMFBoxes(destBoxes.size());  // These are the ones
						  // we want to fillpatch.
    Vector< Vector< Vector< Vector<FillBoxId> > > > fillBoxId;
    Vector< Vector< Vector< Vector<BoxArray> > > >  fillBoxIdBAs;
			          // [grid][level][fillablesubbox][oldnew]
			          // oldnew not used here
    Vector< Vector< Vector<Box> > > savedFineBox;  // [grid][level][fillablesubbox]

    fillBoxId.resize(destBoxes.size());
    fillBoxIdBAs.resize(destBoxes.size());
    savedFineBox.resize(destBoxes.size());
    for(int iBox(0); iBox < destBoxes.size(); ++iBox) {
      if(destMultiFab.DistributionMap()[iBox] == myProc) {
	localMFBoxes.set(iBox, destBoxes[iBox]);
        fillBoxId[iBox].resize(finestFillLevel + 1);
        fillBoxIdBAs[iBox].resize(finestFillLevel + 1);
        savedFineBox[iBox].resize(finestFillLevel + 1);
      }
    }

    IndexType boxType(destBoxes.ixType());
    BoxList unfilledBoxesOnThisLevel(boxType);
    BoxList unfillableBoxesOnThisLevel(boxType);
    // Do this for all local fab boxes.
    for(int ibox(0); ibox < localMFBoxes.size(); ++ibox) {
      if(destMultiFab.DistributionMap()[ibox] != myProc) {
	continue;
      }
        unfilledBoxesOnThisLevel.clear();
        BL_ASSERT(unfilledBoxesOnThisLevel.ixType() == boxType);
        BL_ASSERT(unfilledBoxesOnThisLevel.ixType() == localMFBoxes[ibox].ixType());
        unfilledBoxesOnThisLevel.push_back(localMFBoxes[ibox]);
        // Find the boxes that can be filled on each level--these are all
        // defined at their level of refinement.
        bool needsFilling(true);
        for(currentLevel = finestFillLevel; currentLevel >= 0 && needsFilling;
            --currentLevel)
        {
            unfillableBoxesOnThisLevel.clear();
            const Box &currentPDomain = probDomain[currentLevel];

	    int ufbLength(unfilledBoxesOnThisLevel.size());
            fillBoxId[ibox][currentLevel].resize(ufbLength);
            fillBoxIdBAs[ibox][currentLevel].resize(ufbLength);
            savedFineBox[ibox][currentLevel].resize(ufbLength);

            int currentBLI(0);
            for(BoxList::iterator bli = unfilledBoxesOnThisLevel.begin();
	        bli != unfilledBoxesOnThisLevel.end(); ++bli)
	    {
                BL_ASSERT(bli->ok());
                Box coarseDestBox(*bli);
                Box fineTruncDestBox(coarseDestBox & currentPDomain);
                if(fineTruncDestBox.ok()) {
                  fineTruncDestBox.refine(cumulativeRefRatios[currentLevel]);
                  Box tempCoarseBox;
                  if(currentLevel == finestFillLevel) {
                    tempCoarseBox = fineTruncDestBox;
                  } else {
                    tempCoarseBox = fineTruncDestBox;
		    // check this vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                    tempCoarseBox.coarsen(cumulativeRefRatios[currentLevel]);
                  }

                  savedFineBox[ibox][currentLevel][currentBLI] = fineTruncDestBox;
                  BL_ASSERT(localMFBoxes[ibox].intersects(fineTruncDestBox));

                  BoxList tempUnfillableBoxes(boxType);
                  fillBoxId[ibox][currentLevel][currentBLI].resize(1);
                  fillBoxIdBAs[ibox][currentLevel][currentBLI].resize(1);

                  fillBoxId[ibox][currentLevel][currentBLI][0] = 
		      multiFabCopyDesc.AddBox(stateDataMFId[currentLevel],
					      tempCoarseBox, &tempUnfillableBoxes,
					      srcComp, 0, 1);

                  fillBoxIdBAs[ibox][currentLevel][currentBLI][0] =
                      BoxArray(amrex::complementIn(tempCoarseBox,
		               tempUnfillableBoxes));

                  unfillableBoxesOnThisLevel.catenate(tempUnfillableBoxes);
                  ++currentBLI;
                }
            }

            unfilledBoxesOnThisLevel.clear();
            unfilledBoxesOnThisLevel =
                unfillableBoxesOnThisLevel.intersect(currentPDomain);

            if(unfilledBoxesOnThisLevel.isEmpty()) {
              needsFilling = false;
            } else {
              Box coarseLocalMFBox(localMFBoxes[ibox]);
              coarseLocalMFBox.coarsen(cumulativeRefRatios[currentLevel]);
              unfilledBoxesOnThisLevel.intersect(coarseLocalMFBox);
	      // check this vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	      if(currentLevel != 0) {
                unfilledBoxesOnThisLevel.coarsen(refRatio[currentLevel - 1]);
	      }

              if(currentLevel == 0) {
                BoxList unfilledInside =
                        unfilledBoxesOnThisLevel.intersect(currentPDomain);
                if( ! unfilledInside.isEmpty()) {
                  unfilledInside.intersect(coarseLocalMFBox);
                  BL_ASSERT(unfilledInside.isEmpty());
                }
              }
            }
        }
    }

    multiFabCopyDesc.CollectData();


    for(int currentIndex = 0; currentIndex < destBoxes.size(); ++currentIndex) {
      if(destMultiFab.DistributionMap()[currentIndex] != myProc) {
	continue;
      }
      for(int currLevel(0); currLevel <= finestFillLevel; ++currLevel) {
        for(int currentBox(0);
            currentBox < fillBoxId[currentIndex][currLevel].size();
            ++currentBox)
        {
            Box tempCoarseBox(
		       fillBoxId[currentIndex][currLevel][currentBox][0].box());
            FArrayBox tempCoarseDestFab(tempCoarseBox, 1);
            tempCoarseDestFab.setVal<RunOn::Host>(1.e30);
            multiFabCopyDesc.FillFab(stateDataMFId[currLevel],
			  fillBoxId[currentIndex][currLevel][currentBox][0],
			  tempCoarseDestFab);

            Box intersectDestBox(savedFineBox[currentIndex][currLevel][currentBox]);
            intersectDestBox &= destMultiFab[currentIndex].box();

            const BoxArray &filledBoxes =
                fillBoxIdBAs[currentIndex][currLevel][currentBox][0];
            BoxArray fboxes(filledBoxes);
            FArrayBox *copyFromThisFab;
            const BoxArray *copyFromTheseBoxes;
            FArrayBox tempCurrentFillPatchedFab;

            if(intersectDestBox.ok()) {
              if(currLevel != finestFillLevel) {
                fboxes.refine(cumulativeRefRatios[currLevel]);
                // Interpolate up to fine patch.
                tempCurrentFillPatchedFab.resize(intersectDestBox, nFillComps);
                tempCurrentFillPatchedFab.setVal<RunOn::Host>(1.e30);
		BL_ASSERT(intersectDestBox.ok());
		BL_ASSERT(tempCoarseDestFab.box().ok());
		PcInterp(tempCurrentFillPatchedFab,
			 tempCoarseDestFab, intersectDestBox,
			 cumulativeRefRatios[currLevel]);
                copyFromThisFab = &tempCurrentFillPatchedFab;
                copyFromTheseBoxes = &fboxes;
              } else {
                copyFromThisFab = &tempCoarseDestFab;
                copyFromTheseBoxes = &filledBoxes;
              }
              for(int iFillBox(0); iFillBox < copyFromTheseBoxes->size();
                  ++iFillBox)
              {
                Box srcdestBox((*copyFromTheseBoxes)[iFillBox]);
                srcdestBox &= destMultiFab[currentIndex].box();
                srcdestBox &= intersectDestBox;
                if(srcdestBox.ok()) {
                  destMultiFab[currentIndex].copy<RunOn::Host>(*copyFromThisFab,
                                                  srcdestBox, 0, srcdestBox,
                                                  destComp, nFillComps);
                }
              }
            }
        }
      }  // end for(currentLevel...)
    }  // end for(currentIndex...)

  }  // end for(currentFillIndex...)
}


// ---------------------------------------------------------------
void AmrData::FillVar(Vector<FArrayBox *> &destFabs, const Vector<Box> &destBoxes,
		      int finestFillLevel, const string &varname, int procWithFabs)
{

//
// This function fills dest only on procWithFabs.  All other dest
// pointers (on other processors) should be NULL.  destBox
// on all processors must be defined.
//

   BL_ASSERT(finestFillLevel >= 0 && finestFillLevel <= finestLevel);
   BL_ASSERT(procWithFabs >= 0 && procWithFabs < ParallelDescriptor::NProcs());
   for(int iIndex(0); iIndex < destBoxes.size(); ++iIndex) {
     BL_ASSERT(probDomain[finestFillLevel].contains(destBoxes[iIndex]));
   }

    int myproc(ParallelDescriptor::MyProc());
    int stateIndex(StateNumber(varname));
    int srcComp(0);
    int destComp(0);
    int numFillComps(1);

    int currentLevel;
    Vector<int> cumulativeRefRatios(finestFillLevel + 1, -1);

    cumulativeRefRatios[finestFillLevel] = 1;
    for(currentLevel = finestFillLevel - 1; currentLevel >= 0; --currentLevel) {
      cumulativeRefRatios[currentLevel] = cumulativeRefRatios[currentLevel + 1] *
                                          refRatio[currentLevel];
    }

    // ensure the required grids are in memory
    for(currentLevel = 0; currentLevel <= finestFillLevel; ++currentLevel) {
      for(int iBox = 0; iBox < destBoxes.size(); ++iBox) {
	Box tempCoarseBox(destBoxes[iBox]);
        if(currentLevel != finestFillLevel) {
          tempCoarseBox.coarsen(cumulativeRefRatios[currentLevel]);
        }
        GetGrids(currentLevel, stateIndex, tempCoarseBox);
      }
    }

    MultiFabCopyDescriptor multiFabCopyDesc;
    Vector<MultiFabId> stateDataMFId(finestFillLevel + 1);
    for(currentLevel = 0; currentLevel <= finestFillLevel; ++currentLevel) {
      stateDataMFId[currentLevel] =
           multiFabCopyDesc.RegisterFabArray(dataGrids[currentLevel][stateIndex]);
    }

    Vector<Box> localMFBoxes;      // These are the ones we want to fillpatch.
    Vector< Vector< Vector< Vector<FillBoxId> > > > fillBoxId;
    Vector< Vector< Vector< Vector<BoxArray> > > >  fillBoxIdBAs;
			          // [grid][level][fillablesubbox][oldnew]
			          // oldnew not used here
    Vector< Vector< Vector<Box> > > savedFineBox;  // [grid][level][fillablesubbox]
    if(myproc == procWithFabs) {
      localMFBoxes = destBoxes;
      fillBoxId.resize(destBoxes.size());
      fillBoxIdBAs.resize(destBoxes.size());
      savedFineBox.resize(destBoxes.size());
      for(int iLocal = 0; iLocal < localMFBoxes.size(); ++iLocal) {
        fillBoxId[iLocal].resize(finestFillLevel + 1);
        fillBoxIdBAs[iLocal].resize(finestFillLevel + 1);
        savedFineBox[iLocal].resize(finestFillLevel + 1);
      }
    }

    IndexType boxType(destBoxes[0].ixType());
    BoxList unfilledBoxesOnThisLevel(boxType);
    BoxList unfillableBoxesOnThisLevel(boxType);
    // Do this for all local fab boxes.
    for(int ibox(0); ibox < localMFBoxes.size(); ++ibox) {
        unfilledBoxesOnThisLevel.clear();
        BL_ASSERT(unfilledBoxesOnThisLevel.ixType() == boxType);
        BL_ASSERT(unfilledBoxesOnThisLevel.ixType() == localMFBoxes[ibox].ixType());
        unfilledBoxesOnThisLevel.push_back(localMFBoxes[ibox]);
        // Find the boxes that can be filled on each level--these are all
        // defined at their level of refinement.
        bool needsFilling(true);
        for(currentLevel = finestFillLevel; currentLevel >= 0 && needsFilling;
            --currentLevel)
        {
            unfillableBoxesOnThisLevel.clear();
            const Box &currentPDomain = probDomain[currentLevel];

	    int ufbLength = unfilledBoxesOnThisLevel.size();
            fillBoxId[ibox][currentLevel].resize(ufbLength);
            fillBoxIdBAs[ibox][currentLevel].resize(ufbLength);
            savedFineBox[ibox][currentLevel].resize(ufbLength);

            int currentBLI(0);
            for(BoxList::iterator bli = unfilledBoxesOnThisLevel.begin();
	        bli != unfilledBoxesOnThisLevel.end(); ++bli)
	    {
                BL_ASSERT(bli->ok());
                Box coarseDestBox(*bli);
                Box fineTruncDestBox(coarseDestBox & currentPDomain);
                if(fineTruncDestBox.ok()) {
                  fineTruncDestBox.refine(cumulativeRefRatios[currentLevel]);
                  Box tempCoarseBox;
                  if(currentLevel == finestFillLevel) {
                    tempCoarseBox = fineTruncDestBox;
                  } else {
                    tempCoarseBox = fineTruncDestBox;
		    // check this vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
                    tempCoarseBox.coarsen(cumulativeRefRatios[currentLevel]);
                  }

                  savedFineBox[ibox][currentLevel][currentBLI] = fineTruncDestBox;
                  BL_ASSERT(localMFBoxes[ibox].intersects(fineTruncDestBox));

                  BoxList tempUnfillableBoxes(boxType);
                  fillBoxId[ibox][currentLevel][currentBLI].resize(1);
                  fillBoxIdBAs[ibox][currentLevel][currentBLI].resize(1);

                  fillBoxId[ibox][currentLevel][currentBLI][0] = 
		      multiFabCopyDesc.AddBox(stateDataMFId[currentLevel],
					      tempCoarseBox, &tempUnfillableBoxes,
					      srcComp, destComp, numFillComps);

                  fillBoxIdBAs[ibox][currentLevel][currentBLI][0] =
                      BoxArray(amrex::complementIn(tempCoarseBox,
		               tempUnfillableBoxes));

                  unfillableBoxesOnThisLevel.catenate(tempUnfillableBoxes);
                  ++currentBLI;
                }
            }

            unfilledBoxesOnThisLevel.clear();
            unfilledBoxesOnThisLevel =
                unfillableBoxesOnThisLevel.intersect(currentPDomain);

            if(unfilledBoxesOnThisLevel.isEmpty()) {
              needsFilling = false;
            } else {
              Box coarseLocalMFBox(localMFBoxes[ibox]);
              coarseLocalMFBox.coarsen(cumulativeRefRatios[currentLevel]);
              unfilledBoxesOnThisLevel.intersect(coarseLocalMFBox);
	      // check this vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv
	      if(currentLevel != 0) {
                unfilledBoxesOnThisLevel.coarsen(refRatio[currentLevel - 1]);
	      }

              if(currentLevel == 0) {
                BoxList unfilledInside =
                        unfilledBoxesOnThisLevel.intersect(currentPDomain);
                if( ! unfilledInside.isEmpty()) {
                  unfilledInside.intersect(coarseLocalMFBox);
                  BL_ASSERT(unfilledInside.isEmpty());
                }
              }
            }
        }
    }

    multiFabCopyDesc.CollectData();


    for(int currentIndex = 0; currentIndex < destBoxes.size(); ++currentIndex) {
      for(int currLevel = 0; currLevel <= finestFillLevel; ++currLevel) {
	if(myproc != procWithFabs) {
	  break;
	}
        for(int currentBox(0);
            currentBox < fillBoxId[currentIndex][currLevel].size();
            ++currentBox)
        {
            Box tempCoarseBox(
		       fillBoxId[currentIndex][currLevel][currentBox][0].box());
            FArrayBox tempCoarseDestFab(tempCoarseBox, numFillComps);
            tempCoarseDestFab.setVal<RunOn::Host>(1.e30);
            multiFabCopyDesc.FillFab(stateDataMFId[currLevel],
			  fillBoxId[currentIndex][currLevel][currentBox][0],
			  tempCoarseDestFab);

            Box intersectDestBox(savedFineBox[currentIndex][currLevel][currentBox]);
            intersectDestBox &= destFabs[currentIndex]->box();

            const BoxArray &filledBoxes =
                fillBoxIdBAs[currentIndex][currLevel][currentBox][0];
            BoxArray fboxes(filledBoxes);
            FArrayBox *copyFromThisFab;
            const BoxArray *copyFromTheseBoxes;
            FArrayBox tempCurrentFillPatchedFab;

            if(intersectDestBox.ok()) {
              if(currLevel != finestFillLevel) {
                    fboxes.refine(cumulativeRefRatios[currLevel]);
                    // Interpolate up to fine patch.
                    tempCurrentFillPatchedFab.resize(intersectDestBox, numFillComps);
                    tempCurrentFillPatchedFab.setVal<RunOn::Host>(1.e30);
		    BL_ASSERT(intersectDestBox.ok());
		    BL_ASSERT( tempCoarseDestFab.box().ok());
		    PcInterp(tempCurrentFillPatchedFab,
			     tempCoarseDestFab,
			     intersectDestBox,
			     cumulativeRefRatios[currLevel]);
                    copyFromThisFab = &tempCurrentFillPatchedFab;
                    copyFromTheseBoxes = &fboxes;
              } else {
                    copyFromThisFab = &tempCoarseDestFab;
                    copyFromTheseBoxes = &filledBoxes;
              }
              for(int iFillBox(0); iFillBox < copyFromTheseBoxes->size();
                  ++iFillBox)
              {
                    Box srcdestBox((*copyFromTheseBoxes)[iFillBox]);
                    srcdestBox &= destFabs[currentIndex]->box();
                    srcdestBox &= intersectDestBox;
                    if(srcdestBox.ok()) {
                        destFabs[currentIndex]->copy<RunOn::Host>(*copyFromThisFab,
                                                   srcdestBox, 0, srcdestBox,
                                                   destComp, numFillComps);
                    }
              }
            }
        }
      }  // end for(currentLevel...)
    }  // end for(currentIndex...)
}    // end FillVar for a fab on a single processor


// ---------------------------------------------------------------
void AmrData::FillInterior(FArrayBox &/*dest*/, int /*level*/, const Box &/*subbox*/) {
   amrex::Abort("Error:  should not be in AmrData::FillInterior");
}


// ---------------------------------------------------------------
int AmrData::NumDeriveFunc() const {
  return (plotVars.size());
}


// ---------------------------------------------------------------
// return true if the given name is the name of a plot variable
// that can be derived from what is known.
// ---------------------------------------------------------------
bool AmrData::CanDerive(const string &name) const {
  for(int i(0); i < plotVars.size(); ++i) {
    if(plotVars[i] == name) {
      return true;
    }
  }
  return false;
}


// ---------------------------------------------------------------
bool AmrData::CanDerive(const Vector<string> &names) const {
  int nFound(0);
  for(int n(0); n < names.size(); ++n) {
    for(int i(0); i < plotVars.size(); ++i) {
      if(plotVars[i] == names[n]) {
        ++nFound;
        continue;
      }
    }
  }
  return(nFound == names.size());
}


// ---------------------------------------------------------------
// output the list of variables that can be derived
// ---------------------------------------------------------------
void AmrData::ListDeriveFunc(std::ostream &os) const {
   for(int i(0); i < plotVars.size(); ++i) {
     os << plotVars[i] << endl;
   }
}


// ---------------------------------------------------------------
int AmrData::NIntersectingGrids(int level, const Box &b) const {
  BL_ASSERT(level >=0 && level <= finestLevel);
  BL_ASSERT(b.ok());

  int nGrids(0);
  if(fileType == Amrvis::FAB || (fileType == Amrvis::MULTIFAB && level == 0)) {
    nGrids = 1;
  } else {
    const BoxArray &visMFBA = visMF[level][0]->boxArray();
    for(int boxIndex = 0; boxIndex < visMFBA.size(); ++boxIndex) {
      if(b.intersects(visMFBA[boxIndex])) {
        ++nGrids;
      }
    }
  }
  return nGrids;
}


// ---------------------------------------------------------------
int AmrData::FinestContainingLevel(const Box &b, int startLevel) const {
  BL_ASSERT(startLevel >= 0 && startLevel <= finestLevel);
  BL_ASSERT(b.ok());

  if(fileType == Amrvis::FAB) {
    return 0;
  } else {
    Box levelBox(b);
    for(int level = startLevel; level > 0; --level) {
      const BoxArray &visMFBA = visMF[level][0]->boxArray();
      if(visMFBA.contains(levelBox)) {
        return level;
      }
      levelBox.coarsen(refRatio[level - 1]);
    }
  }
  return 0;
}


// ---------------------------------------------------------------
int AmrData::FinestIntersectingLevel(const Box &b, int startLevel) const {
  BL_ASSERT(startLevel >= 0 && startLevel <= finestLevel);
  BL_ASSERT(b.ok());

  if(fileType == Amrvis::FAB) {
    return 0;
  } else {
    Box levelBox(b);
    for(int level(startLevel); level > 0; --level) {
      const BoxArray &visMFBA = visMF[level][0]->boxArray();

      for(int iBox(0); iBox < visMFBA.size(); ++iBox) {
        if(visMFBA[iBox].intersects(levelBox)) {
          return level;
	}
      }

      levelBox.coarsen(refRatio[level - 1]);
    }
  }
  return 0;
}


// ---------------------------------------------------------------
MultiFab &AmrData::GetGrids(int level, int componentIndex) {
  for(MFIter mfi(*dataGrids[level][componentIndex]); mfi.isValid(); ++mfi) {
    DefineFab(level, componentIndex, mfi.index());
  }
  return *dataGrids[level][componentIndex];
}


// ---------------------------------------------------------------
MultiFab &AmrData::GetGrids(int level, int componentIndex, const Box &onBox) {
  if(fileType == Amrvis::FAB || (fileType == Amrvis::MULTIFAB && level == 0)) {
    // do nothing
  } else {
    int whichVisMF(compIndexToVisMFMap[componentIndex]);
    for(MFIter mfi(*dataGrids[level][componentIndex]);
        mfi.isValid(); ++mfi)
    {
      if(onBox.intersects(visMF[level][whichVisMF]->boxArray()[mfi.index()])) {
        DefineFab(level, componentIndex, mfi.index());
      }
    }
  }
  return *dataGrids[level][componentIndex];
}


// ---------------------------------------------------------------
bool AmrData::DefineFab(int level, int componentIndex, int fabIndex) {

  if( ! dataGridsDefined[level][componentIndex][fabIndex]) {
    int whichVisMF(compIndexToVisMFMap[componentIndex]);
    int whichVisMFComponent(compIndexToVisMFComponentMap[componentIndex]);
    dataGrids[level][componentIndex]->setFab(fabIndex,
                visMF[level][whichVisMF]->readFAB(fabIndex, whichVisMFComponent));
    dataGridsDefined[level][componentIndex][fabIndex] = true;
  }
  return true;
}


// ---------------------------------------------------------------
void AmrData::FlushGrids() {
  for (int componentIndex(0); componentIndex < nComp; ++componentIndex) {
    FlushGrids(componentIndex);
  }
}


// ---------------------------------------------------------------
void AmrData::FlushGrids(int componentIndex) {

  MFInfo Fab_noallocate;
  Fab_noallocate.SetAlloc(false);

  BL_ASSERT(componentIndex < nComp);
  for(int lev(0); lev <= finestLevel; ++lev) {
    if(dataGrids.size() > lev
       && dataGrids[lev].size() > componentIndex
       && dataGrids[lev][componentIndex])
    {
      BoxArray ba = dataGrids[lev][componentIndex]->boxArray();
      DistributionMapping dm = dataGrids[lev][componentIndex]->DistributionMap();
      int flushNGrow = dataGrids[lev][componentIndex]->nGrow();
      delete dataGrids[lev][componentIndex];
      dataGrids[lev][componentIndex] = new MultiFab(ba, dm, 1, flushNGrow, Fab_noallocate);
      for(MFIter mfi(*dataGrids[lev][componentIndex]); mfi.isValid(); ++mfi) {
         dataGridsDefined[lev][componentIndex][mfi.index()] = false;
      }
    }
  }
}


// ---------------------------------------------------------------
bool AmrData::MinMax(const Box &onBox, const string &derived, int level,
		     Real &dataMin, Real &dataMax)
{
  BL_ASSERT(level >= 0 && level <= finestLevel);
  BL_ASSERT(onBox.ok());

  bool valid(false);  // does onBox intersect any grids (are minmax valid)
  Real minVal, maxVal;
  dataMin =  std::numeric_limits<Real>::max();
  dataMax = -std::numeric_limits<Real>::max();
  Box overlap;

  //  our strategy here is to use the VisMF min and maxes if possible
  //  first, test if onBox completely contains each multifab box
  //  if so, use VisMF min and max
  //  if not, test if VisMF min and max are within dataMin and dataMax
  //  if so, use VisMF min and max

  int compIndex(StateNumber(derived));

  if(fileType == Amrvis::FAB || (fileType == Amrvis::MULTIFAB && level == 0)) {
    for(MFIter gpli(*dataGrids[level][compIndex]); gpli.isValid(); ++gpli) {
      if(onBox.intersects(dataGrids[level][compIndex]->boxArray()[gpli.index()])) {
          valid = true;
          overlap = onBox;
          overlap &= gpli.validbox();
          minVal = (*dataGrids[level][compIndex])[gpli].min<RunOn::Host>(overlap, 0);
          maxVal = (*dataGrids[level][compIndex])[gpli].max<RunOn::Host>(overlap, 0);

          dataMin = std::min(dataMin, minVal);
          dataMax = std::max(dataMax, maxVal);
      }
    }

  } else if(bCartGrid && (compIndex != StateNumber("vfrac")) && bShowBody) {
#if (BL_SPACEDIM == 1)
    amrex::Abort("AmrData::MinMax:  should not be here for 1d.");
#else
    int iCount(0), iCountAllBody(0);
    int iCountMixedSkipped(0), iCountMixedFort(0);
    int cCount(0), cCountAllBody(0), cCountAllFluid(0), cCountMixed(0);
    int cCountMixedSkipped(0), cCountMixedFort(0);
    int allCount(0), outsideCount(0);
    for(MFIter gpli(*dataGrids[level][compIndex]); gpli.isValid(); ++gpli) {
      ++allCount;
      int gdx(gpli.index());
      int whichVisMF(compIndexToVisMFMap[compIndex]);
      int whichVisMFComponent(compIndexToVisMFComponentMap[compIndex]);
      Real visMFMin(visMF[level][whichVisMF]->min(gdx, whichVisMFComponent));
      Real visMFMax(visMF[level][whichVisMF]->max(gdx, whichVisMFComponent));
      int vfIndex(StateNumber("vfrac"));
      if(onBox.contains(gpli.validbox())) {
        ++cCount;
        int vfWhichVisMF(compIndexToVisMFMap[vfIndex]);
        int vfWhichVisMFComponent(compIndexToVisMFComponentMap[vfIndex]);
        Real vfVisMFMin(visMF[level][vfWhichVisMF]->min(gdx, vfWhichVisMFComponent));
        Real vfVisMFMax(visMF[level][vfWhichVisMF]->max(gdx, vfWhichVisMFComponent));
        if(vfVisMFMin > (1.0 - vfEps[level])) {  // no cg body in this grid
	  ++cCountAllFluid;
          dataMin = std::min(dataMin, visMFMin);
          dataMax = std::max(dataMax, visMFMax);
          valid = true;
	} else if(vfVisMFMax >= vfEps[level] ) {
	  ++cCountMixed;
          if(visMFMin < dataMin || visMFMax > dataMax) {  // do it the hard way
            DefineFab(level, compIndex, gdx);
            DefineFab(level, vfIndex, gdx);
            Real *ddat = (*dataGrids[level][compIndex])[gpli].dataPtr();
            Real *vdat = (*dataGrids[level][vfIndex])[gpli].dataPtr();
            const int *dlo = (*dataGrids[level][compIndex])[gpli].loVect();
            const int *dhi = (*dataGrids[level][compIndex])[gpli].hiVect();

            overlap = onBox;
            overlap &= gpli.validbox();
            Real vfMaxVal = (*dataGrids[level][vfIndex])[gpli].max<RunOn::Host>(overlap, 0);
            if(vfMaxVal >= vfEps[level]) {
	      ++cCountMixedFort;
              valid = true;

              FORT_CARTGRIDMINMAX(ddat, AMREX_ARLIM(dlo), AMREX_ARLIM(dhi), vdat, vfEps[level],
                                  minVal, maxVal);
              dataMin = std::min(dataMin, minVal);
              dataMax = std::max(dataMax, maxVal);
            }
	  } else {
	    ++cCountMixedSkipped;
	  }
	} else {  // all body
	  ++cCountAllBody;
	}
      } else if(onBox.intersects(visMF[level][whichVisMF]->boxArray()[gdx])) {
	++iCount;
        if(visMFMin < dataMin || visMFMax > dataMax) {  // do it the hard way
          DefineFab(level, compIndex, gdx);
          DefineFab(level, vfIndex, gdx);
          Real *ddat = (*dataGrids[level][compIndex])[gpli].dataPtr();
          Real *vdat = (*dataGrids[level][vfIndex])[gpli].dataPtr();
          const int *dlo = (*dataGrids[level][compIndex])[gpli].loVect();
          const int *dhi = (*dataGrids[level][compIndex])[gpli].hiVect();

          overlap = onBox;
          overlap &= gpli.validbox();
          Real vfMaxVal = (*dataGrids[level][vfIndex])[gpli].max<RunOn::Host>(overlap, 0);
          if(vfMaxVal >= vfEps[level]) {
	    ++iCountMixedFort;
            valid = true;

            FORT_CARTGRIDMINMAX(ddat, AMREX_ARLIM(dlo), AMREX_ARLIM(dhi), vdat, vfEps[level],
                                minVal, maxVal);
            dataMin = std::min(dataMin, minVal);
            dataMax = std::max(dataMax, maxVal);
          } else {
	    ++iCountAllBody;
	  }
	} else {
	  ++iCountMixedSkipped;
	}
      } else {
        ++outsideCount;
      }
    }
    if(verbose) {
      cout << "*** Level:  " << level << endl;
      SHOWVAL(allCount);
      SHOWVAL(outsideCount);
      cout << "--------------" << endl;

      SHOWVAL(cCount);
      SHOWVAL(cCountAllBody);
      SHOWVAL(cCountAllFluid);
      SHOWVAL(cCountMixed);
      SHOWVAL(cCountMixedSkipped);
      SHOWVAL(cCountMixedFort);
      cout << "--------------" << endl;
      if(iCount > 0) {
        SHOWVAL(iCount);
        SHOWVAL(iCountAllBody);
        SHOWVAL(iCountMixedSkipped);
        SHOWVAL(iCountMixedFort);
        cout << "--------------" << endl;
      }
      cout << endl << endl;
    }
#endif

  } else {
    for(MFIter gpli(*dataGrids[level][compIndex]); gpli.isValid(); ++gpli) {
      int whichVisMF(compIndexToVisMFMap[compIndex]);
      int whichVisMFComponent(compIndexToVisMFComponentMap[compIndex]);
      Real visMFMin(visMF[level][whichVisMF]->min(gpli.index(),
		    whichVisMFComponent));
      Real visMFMax(visMF[level][whichVisMF]->max(gpli.index(),
		    whichVisMFComponent));
      if(onBox.contains(gpli.validbox())) {
	  dataMin = std::min(dataMin, visMFMin);
	  dataMax = std::max(dataMax, visMFMax);
	  valid = true;
      } else if(onBox.intersects(visMF[level][whichVisMF]->
				 boxArray()[gpli.index()]))
      {
        if(visMFMin < dataMin || visMFMax > dataMax) {  // do it the hard way
	  DefineFab(level, compIndex, gpli.index());
          valid = true;
          overlap = onBox;
          overlap &= gpli.validbox();
          minVal = (*dataGrids[level][compIndex])[gpli].min<RunOn::Host>(overlap, 0);
          maxVal = (*dataGrids[level][compIndex])[gpli].max<RunOn::Host>(overlap, 0);

          dataMin = std::min(dataMin, minVal);
          dataMax = std::max(dataMax, maxVal);
        }  // end if(visMFMin...)
      }
    }
  }

  ParallelDescriptor::ReduceRealMin(dataMin);
  ParallelDescriptor::ReduceRealMax(dataMax);

  return valid;
}  // end MinMax


// ---------------------------------------------------------------
int AmrData::StateNumber(const string &statename) const {
  for(int ivar(0); ivar < plotVars.size(); ++ivar) {
    if(statename == plotVars[ivar]) {
      return ivar;
    }
  }
  if(ParallelDescriptor::IOProcessor()) {
    cerr << "Error:  bad state name:  " << statename << std::endl;
  }
  amrex::Abort("bad state name in AmrData::StateNumber()");
  return(-1);
}


// ---------------------------------------------------------------
void AmrData::Interp(FArrayBox &fine, FArrayBox &crse,
                     const Box &fine_box, int lrat)
{
#if (BL_SPACEDIM == 1)
    amrex::Abort("AmrData::MinMax:  should not be here for 1d.");
#else
   BL_ASSERT(fine.box().contains(fine_box));
   Box crse_bx(amrex::coarsen(fine_box,lrat));
   Box fslope_bx(amrex::refine(crse_bx,lrat));
   Box cslope_bx(crse_bx);
   cslope_bx.grow(1);
   BL_ASSERT(crse.box() == cslope_bx);

   // alloc temp space for coarse grid slopes
   Long cLen = cslope_bx.numPts();
   Real *cslope = new Real[BL_SPACEDIM*cLen];
   Long loslp    = cslope_bx.index(crse_bx.smallEnd());
   Long hislp    = cslope_bx.index(crse_bx.bigEnd());
   Long cslope_vol = cslope_bx.numPts();
   Long clo = 1 - loslp;
   Long chi = clo + cslope_vol - 1;
   cLen = hislp - loslp + 1;

   // alloc temp space for one strip of fine grid slopes
   int dir;
   int fLen = fslope_bx.longside(dir);
   Real *fdat   = new Real[(BL_SPACEDIM+2)*fLen];
   Real *foff   = fdat + fLen;
   Real *fslope = foff + fLen;


   // alloc tmp space for slope calc and to allow for vectorization
   const int *fblo = fine_box.loVect();
   const int *fbhi = fine_box.hiVect();
   const int *cblo = crse_bx.loVect();
   const int *cbhi = crse_bx.hiVect();
   const int *fslo = fslope_bx.loVect();
   const int *fshi = fslope_bx.hiVect();

   FORT_CINTERP(fine.dataPtr(0),AMREX_ARLIM(fine.loVect()),AMREX_ARLIM(fine.hiVect()),
               fblo,fbhi,fine.nComp(),lrat,
               crse.dataPtr(0),clo,chi,cblo,cbhi,fslo,fshi,
               cslope,cLen,fslope,fdat,fLen,foff);

   delete [] fdat;
   delete [] cslope;
#endif
}


// ---------------------------------------------------------------
void AmrData::PcInterp(FArrayBox &fine, const FArrayBox &crse,
                       const Box &subbox, int lrat)
{
   BL_ASSERT(fine.box().contains(subbox));
   BL_ASSERT(fine.nComp() == crse.nComp());
   Box cfine(crse.box());
   cfine.refine(lrat);
   Box fine_ovlp(subbox);
   fine_ovlp &= cfine;
   if(fine_ovlp.ok()) {
      const int *fblo = fine_ovlp.smallEnd().getVect();
      const int *fbhi = fine_ovlp.bigEnd().getVect();
      Box crse_ovlp(fine_ovlp);
      crse_ovlp.coarsen(lrat);
      const int *cblo = crse_ovlp.smallEnd().getVect();
      const int *cbhi = crse_ovlp.bigEnd().getVect();
      Box fine_temp(crse_ovlp);
      fine_temp.refine(lrat);
      int tlo = fine_temp.smallEnd()[0];
      int thi = fine_temp.bigEnd()[0];
      int inextra(0);
      if(fine_temp.ixType().test(0) == true) {  // node type
        inextra = 1;
      }
      Real *tempSpace = new Real[thi-tlo+1+inextra];
      FORT_PCINTERP(fine.dataPtr(0),AMREX_ARLIM(fine.loVect()),AMREX_ARLIM(fine.hiVect()),
                   fblo,fbhi, lrat,fine.nComp(),
                   crse.dataPtr(),AMREX_ARLIM(crse.loVect()),AMREX_ARLIM(crse.hiVect()),
                   cblo,cbhi, tempSpace,tlo,thi);

      delete [] tempSpace;
   }
}


// ---------------------------------------------------------------
FArrayBox *AmrData::ReadGrid(std::istream &is, int numVar) {
   Long i, gstep;
   Real timeIn;
   static int gridCount(0);
   Box gbox;
   int glev;

   int gid(gridCount);
   ++gridCount;

   is >> gbox >> glev;
   VSHOWVAL(verbose, gbox)
   VSHOWVAL(verbose, glev)

   is >> gstep >> timeIn;
   VSHOWVAL(verbose, gstep)
   VSHOWVAL(verbose, timeIn)

   for(i = 0; i < BL_SPACEDIM; ++i) {
     Real xlo, xhi;
     is >> xlo >> xhi;  // unused
     if(verbose) {
       cout << "xlo xhi [" << i << "] = " << xlo << "  " << xhi << endl;
     }
   }
   while (is.get() != '\n') {
     ;  // do nothing
   }

   FArrayBox *fabPtr = new FArrayBox(gbox, numVar);
   int whileTrap(0);
   int ivar(0);
   //  optimize this for numVar == newdat.nComp()
   while(ivar < numVar) {
     //FArrayBox tempfab(is);
     FArrayBox tempfab;
     tempfab.readFrom(is);
     fabPtr->copy<RunOn::Host>(tempfab, 0, ivar, tempfab.nComp());
     ivar += tempfab.nComp();
     if(++whileTrap > 256) {   // an arbitrarily large number
       cerr << "Error in GridPlot:  whileTrap caught loop." << endl;
       exit(-4);
     }
   }

   if(verbose) {
     cout << "Constructing Grid, lev = " << glev << "  id = " << gid;
     cout << " box = " << gbox << endl;
   }
  return fabPtr;
}
// ---------------------------------------------------------------
// ---------------------------------------------------------------

}
