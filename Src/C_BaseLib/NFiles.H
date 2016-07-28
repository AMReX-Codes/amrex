
#ifndef BL_NFILES_H
#define BL_NFILES_H

#include <string>
#include <fstream>

#include <ParallelDescriptor.H>

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

    // ---- the nfiles will be named filePrefix00000 through
    // ---- filePrefixnoutfile-1
    // ---- the range [1, nProcs-1] for noutfiles is enforced
    NFilesIter(int noutfiles, const std::string &filePrefix);
    ~NFilesIter();
    bool ReadyToWrite();
    NFilesIter &operator++();

    std::ofstream &Stream() { return fileStream; }

    // ---- get the current Stream()'s seek position
    std::streampos SeekPos();

    static int FileNumber(int nOutFiles, int whichProc)
    {
      return(whichProc % nOutFiles);
    }

    static std::string FileName(int nOutFiles,
                                const std::string &filePrefix,
				int whichProc)
    {
      return BoxLib::Concatenate(filePrefix, FileNumber(nOutFiles, whichProc), 5);
    }


  private:

    int myProc;
    int nProcs;
    int nOutFiles;
    int nSets;
    int mySet;
    std::string fullFileName;
    std::ofstream fileStream;
    bool finishedWriting;

    NFilesIter();  // ---- disallow
};

#endif  /* BL_NFILES_H */
