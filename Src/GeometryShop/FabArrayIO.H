/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _FabArrayIO__H_
#define _FabArrayIO__H_


#include <cmath>
#include <cstdlib>

#include "AMReX_REAL.H"
#include "AMReX_LoHiSide.H"
#include "AMReX_RealVect.H"
#include "AMReX_Box.H"
#include "AMReX_IntVect.H"
#include "AMReX_EBISBox.H"
#include "AMReX_VolIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_IrregNode.H"


namespace amrex
{
  template <class T >
  long
  FabArrayWrite (const FabArray<T>&    mf,
                 const std::string& mf_name)
{
    BL_PROFILE("VisMF::Write_FabArray");
    BL_ASSERT(mf_name[mf_name.length() - 1] != '/');

    // ---- add stream retry
    // ---- add stream buffer (to nfiles)

    int coordinatorProc(ParallelDescriptor::IOProcessorNumber());
    long bytesWritten(0);
    bool calcMinMax(false);
    VisMF::Header hdr(mf, how, currentVersion, calcMinMax);

    std::string filePrefix(mf_name + FabFileSuffix);

    NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);

    bool oldHeader(currentVersion == VisMF::Header::Version_v1);

      if(useDynamicSetSelection) {
        nfi.SetDynamic();
      }
      for( ; nfi.ReadyToWrite(); ++nfi) {
	  // ---- find the total number of bytes including fab headers if needed
          const FABio &fio = FArrayBox::getFABio();
          int whichRDBytes(whichRD->numBytes()), nFABs(0);
          long writeDataItems(0), writeDataSize(0);
          for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	    const FArrayBox &fab = mf[mfi];
	    if(oldHeader) {
	      std::stringstream hss;
	      fio.write_header(hss, fab, fab.nComp());
	      bytesWritten += hss.tellp();
	    }
	    bytesWritten += fab.box().numPts() * mf.nComp() * whichRDBytes;
	    ++nFABs;
	  }
	  char *allFabData(nullptr);
	  bool canCombineFABs(false);
	  if((nFABs > 1 || doConvert) && VisMF::useSingleWrite) {
	    allFabData = new(std::nothrow) char[bytesWritten];
	  }    // ---- else { no need to make a copy for one fab }
	  if(allFabData == nullptr) {
	    canCombineFABs = false;
	  } else {
	    canCombineFABs = true;
	  }

	  if(canCombineFABs) {
            long writePosition(0);
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      char *afPtr = allFabData + writePosition;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
	        memcpy(afPtr, hss.str().c_str(), hLength);  // ---- the fab header
	      }
	      if(doConvert) {
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (afPtr + hLength),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
	      } else {    // ---- copy from the fab
	        memcpy(afPtr + hLength, fab.dataPtr(), writeDataSize);
	      }
              writePosition += hLength + writeDataSize;
            }
            nfi.Stream().write(allFabData, bytesWritten);
            nfi.Stream().flush();
	    delete [] allFabData;

	  } else {    // ---- write fabs individually
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
                nfi.Stream().write(hss.str().c_str(), hLength);    // ---- the fab header
                nfi.Stream().flush();
	      }
	      if(doConvert) {
	        char *cDataPtr = new char[writeDataSize];
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (cDataPtr),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
                nfi.Stream().write(cDataPtr, writeDataSize);
                nfi.Stream().flush();
	        delete [] cDataPtr;
	      } else {    // ---- copy from the fab
                nfi.Stream().write((char *) fab.dataPtr(), writeDataSize);
                nfi.Stream().flush();
	      }
            }
	  }
      }


    if(useDynamicSetSelection) {
      coordinatorProc = nfi.CoordinatorProc();
    }

    if(currentVersion == VisMF::Header::Version_v1 ||
       currentVersion == VisMF::Header::NoFabHeaderMinMax_v1)
    {
      hdr.CalculateMinMax(mf, coordinatorProc);
    }

    VisMF::FindOffsets(mf, filePrefix, hdr, groupSets, currentVersion,
		       useDynamicSetSelection, nfi);

    bytesWritten += VisMF::WriteHeader(mf_name, hdr, coordinatorProc);

    delete whichRD;

    return bytesWritten;
}
    // ---- add stream buffer (to nfiles)
    RealDescriptor *whichRD;
    if(FArrayBox::getFormat() == FABio::FAB_NATIVE) {
      whichRD = FPC::NativeRealDescriptor().clone();
    } else if(FArrayBox::getFormat() == FABio::FAB_NATIVE_32) {
      whichRD = FPC::Native32RealDescriptor().clone();
    } else if(FArrayBox::getFormat() == FABio::FAB_IEEE_32) {
      whichRD = FPC::Ieee32NormalRealDescriptor().clone();
    }
    bool doConvert(*whichRD != FPC::NativeRealDescriptor());

    if(set_ghost) {
        FabArray<FArrayBox>* the_mf = const_cast<FabArray<FArrayBox>*>(&mf);

        for(MFIter mfi(*the_mf); mfi.isValid(); ++mfi) {
            const int idx(mfi.index());

            for(int j(0); j < mf.nComp(); ++j) {
                const Real valMin(mf[mfi].min(mf.box(idx), j));
                const Real valMax(mf[mfi].max(mf.box(idx), j));
                const Real val((valMin + valMax) / 2.0);

                the_mf->get(mfi).setComplement(val, mf.box(idx), j, 1);
            }
        }
    }

    int coordinatorProc(ParallelDescriptor::IOProcessorNumber());
    long bytesWritten(0);
    bool calcMinMax(false);
    VisMF::Header hdr(mf, how, currentVersion, calcMinMax);

    std::string filePrefix(mf_name + FabFileSuffix);

    NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);

    bool oldHeader(currentVersion == VisMF::Header::Version_v1);

      if(useDynamicSetSelection) {
        nfi.SetDynamic();
      }
      for( ; nfi.ReadyToWrite(); ++nfi) {
	  // ---- find the total number of bytes including fab headers if needed
          const FABio &fio = FArrayBox::getFABio();
          int whichRDBytes(whichRD->numBytes()), nFABs(0);
          long writeDataItems(0), writeDataSize(0);
          for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
	    const FArrayBox &fab = mf[mfi];
	    if(oldHeader) {
	      std::stringstream hss;
	      fio.write_header(hss, fab, fab.nComp());
	      bytesWritten += hss.tellp();
	    }
	    bytesWritten += fab.box().numPts() * mf.nComp() * whichRDBytes;
	    ++nFABs;
	  }
	  char *allFabData(nullptr);
	  bool canCombineFABs(false);
	  if((nFABs > 1 || doConvert) && VisMF::useSingleWrite) {
	    allFabData = new(std::nothrow) char[bytesWritten];
	  }    // ---- else { no need to make a copy for one fab }
	  if(allFabData == nullptr) {
	    canCombineFABs = false;
	  } else {
	    canCombineFABs = true;
	  }

	  if(canCombineFABs) {
            long writePosition(0);
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      char *afPtr = allFabData + writePosition;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
	        memcpy(afPtr, hss.str().c_str(), hLength);  // ---- the fab header
	      }
	      if(doConvert) {
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (afPtr + hLength),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
	      } else {    // ---- copy from the fab
	        memcpy(afPtr + hLength, fab.dataPtr(), writeDataSize);
	      }
              writePosition += hLength + writeDataSize;
            }
            nfi.Stream().write(allFabData, bytesWritten);
            nfi.Stream().flush();
	    delete [] allFabData;

	  } else {    // ---- write fabs individually
            for(MFIter mfi(mf); mfi.isValid(); ++mfi) {
              int hLength(0);
              const FArrayBox &fab = mf[mfi];
	      writeDataItems = fab.box().numPts() * mf.nComp();
	      writeDataSize = writeDataItems * whichRDBytes;
	      if(oldHeader) {
	        std::stringstream hss;
	        fio.write_header(hss, fab, fab.nComp());
	        hLength = hss.tellp();
                nfi.Stream().write(hss.str().c_str(), hLength);    // ---- the fab header
                nfi.Stream().flush();
	      }
	      if(doConvert) {
	        char *cDataPtr = new char[writeDataSize];
	        RealDescriptor::convertFromNativeFormat(static_cast<void *> (cDataPtr),
		                                        writeDataItems,
		                                        fab.dataPtr(), *whichRD);
                nfi.Stream().write(cDataPtr, writeDataSize);
                nfi.Stream().flush();
	        delete [] cDataPtr;
	      } else {    // ---- copy from the fab
                nfi.Stream().write((char *) fab.dataPtr(), writeDataSize);
                nfi.Stream().flush();
	      }
            }
	  }
      }


    if(useDynamicSetSelection) {
      coordinatorProc = nfi.CoordinatorProc();
    }

    if(currentVersion == VisMF::Header::Version_v1 ||
       currentVersion == VisMF::Header::NoFabHeaderMinMax_v1)
    {
      hdr.CalculateMinMax(mf, coordinatorProc);
    }

    VisMF::FindOffsets(mf, filePrefix, hdr, groupSets, currentVersion,
		       useDynamicSetSelection, nfi);

    bytesWritten += VisMF::WriteHeader(mf_name, hdr, coordinatorProc);

    delete whichRD;

    return bytesWritten;
}

#endif
