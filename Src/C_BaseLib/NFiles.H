
#ifndef BL_NFILES_H
#define BL_NFILES_H

#include <string>
#include <fstream>

#include <ParallelDescriptor.H>
#include <VisMF.H>

class OrientationIter;

//
// This class encapsulates writing to nfiles with
// static set selection.
//
// an example:
// 
// for(NFilesIter nfi(nOutFiles, filePrefix); nfi.ReadyToWrite(); ++nfi) {
//   // ---- do your writing here
//   nfi.Stream().write((char *) data.dataPtr(), data.size() * sizeof(int));
// }

//

class NFilesIter
{

  public:

    // ---- the nfiles will be named "filePrefix" + "00000"
    // ---- through "filePrefix" + "noutfiles-1"
    // ---- the range [1, nProcs] for noutfiles is enforced
    NFilesIter(int noutfiles, const std::string &filePrefix,
               bool setBuf = false);

    ~NFilesIter();

    bool ReadyToWrite();
    
    NFilesIter &operator++();

    std::ofstream &Stream() { return fileStream; }

    // ---- get the current Stream()'s seek position
    std::streampos SeekPos();

    // ---- this returns the actual number of files used
    // ---- the range [1, nProcs] is enforced
    static int ActualNFiles(int nOutFiles)
    {
      return( std::max(1, std::min(ParallelDescriptor::NProcs(), nOutFiles)) );
    }

    static int FileNumber(int nOutFiles, int whichProc)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());
      return(whichProc % ActualNFiles(nOutFiles));
    }

    static std::string FileName(int nOutFiles,
                                const std::string &filePrefix,
				int whichProc)
    {
      BL_ASSERT(whichProc >= 0 && whichProc < ParallelDescriptor::NProcs());
      return ( BoxLib::Concatenate(filePrefix,
                                   FileNumber(ActualNFiles(nOutFiles),
				              whichProc), 5 ) );
    }


  private:

    int myProc;
    int nProcs;
    int nOutFiles;
    int nSets;
    int mySet;
    std::string fullFileName;
    std::ofstream fileStream;
    VisMF::IO_Buffer io_buffer;
    bool finishedWriting;

    NFilesIter();  // ---- disallow
};

#endif  /* BL_NFILES_H */
