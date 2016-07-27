
#ifndef BL_NFILES_H
#define BL_NFILES_H

#include <string>
#include <fstream>

#include <ParallelDescriptor.H>

class OrientationIter;

//
// Iterate over nfiles
//

class NFilesIter
{

  public:

    NFilesIter(int noutfiles, const std::string &filePrefix);
    ~NFilesIter();
    bool ReadyToWrite();
    NFilesIter &operator++();

    std::ofstream &Stream() { return fileStream; }
    std::streampos SeekPos();


  private:

    int myProc;
    int nProcs;
    int nOutFiles;
    int nSets;
    int mySet;
    int whichSet;
    std::string fullFileName;
    std::ofstream fileStream;
    bool finishedWriting;

    NFilesIter();  // ---- disallow
};

#endif  /* BL_NFILES_H */
