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

#include "AMReX_VisMF.H"
#include "AMReX_EBArith.H"
#include "AMReX_SPMD.H"
#include "AMReX_Print.H"
#include "AMReX_Utility.H"
#include "AMReX_BoxArray.H"


namespace amrex
{
  ///
  class FAIOElement
  {
  public:

    ///
    FAIOElement (int a_fileid, long a_head, int a_boxid);
    FAIOElement ()
      {
      }

    /// needed for sorting
    bool operator<(const FAIOElement& a_fab) const;

    ///
    int linearSize() const;

    ///
    void linearOut(void* buffer ) const;

    ///
    void linearIn(void* buffer );

    string      m_filename;
    int         m_procid;  // used to create the filename
    long        m_head;    // Offset to start of FAB in file.
    int         m_boxid;   // integer into the box array
  };

  ///write a header to   disk in ascii
  std::ostream& 
  operator<< (std::ostream& a_os, FAIOElement& a_elem);

  ///read a header from disk in ascii
  std::istream& 
  operator>> (std::istream& a_is, FAIOElement& a_elem);

  ///
  template<class T>
  class FAIOHeader
  {
  public:
    
    FAIOHeader(const FabArray<T> & a_data)
      {
        m_ba = a_data.boxArray();
        vector<FAIOElement> localElements(a_data.local_size());
        int ielem  = 0;
        int offset = 0;
        int startcomp = 0; int endcomp = a_data.nComp();  


        for(MFIter mfi(a_data); mfi.isValid(); ++mfi) 
        {
          int boxid   = mfi.index();
          //writes one file per proc
          //anything more flexible requires a bunch of communication
          int filenum = procID();
          Box valid   = mfi.validbox();
          int numBytesThisBox = a_data[mfi].nBytes(valid, startcomp, endcomp);
          localElements[ielem] = FAIOElement(filenum, offset, boxid);
          ielem++;
          offset += numBytesThisBox;
        }
        
        vector<vector<FAIOElement> > allElements;
        gather(allElements, localElements, 0);
        if(procID() == 0)
        {
          m_vecElem.resize(0);
          for(int ivec = 0; ivec < allElements.size();  ivec++)
          {
            const vector<FAIOElement>& onevec = allElements[ivec];
            m_vecElem.insert(m_vecElem.end(), onevec.begin(), onevec.end());
          }
        }
        broadcast(m_vecElem, 0);
        std::sort(m_vecElem.begin(), m_vecElem.end()); //this should get them in the same order as the boxarray
      }

    BoxArray                   m_ba;        
    std::vector<FAIOElement>   m_vecElem;   

  private:
    //no weak construction
    FAIOHeader();
  };

  ///write a header to   disk in ascii
  template <class T>
  std::ostream& 
  operator<< (std::ostream& a_os, FAIOHeader<T>& a_hd)
  {
    a_os << a_hd.m_ba         << endl;
    a_os << a_hd.m_vecElem.size() << endl;
    for(int ivec = 0; ivec < a_hd.m_vecElem.size(); ivec++)
    {
      a_os << a_hd.m_vecElem[ivec] << endl;
    }
    return a_os;
  }

  ///read a header from disk in ascii
  template <class T>
  std::istream& 
  operator>> (std::istream& a_is, FAIOHeader<T>& a_hd)
  {
    a_is >> a_hd.m_ba;
    int isize;
    a_is >> isize;
    a_hd.m_vecElem.resize(isize);
    for(int ivec = 0; ivec < a_hd.m_vecElem.size(); ivec++)
    {
      a_is >> a_hd.m_vecElem[ivec];
    }
    return a_is;
  }

  ///class to do I/O of more stuff that is more general than FArrayBox. (reduces everything to streams of bytes)
  /**
     writes one file per proc
     anything more flexible requires a bunch of communication
  */
  template <class T >
    class FabArrayIO
    {
    public:
      ///
      static void
      write(const FabArray<T>&    a_data,
            const std::string&    a_directory_name)
        {
          BL_PROFILE("FabArrayIO::write");
          BL_ASSERT(a_directory_name[a_directory_name.length() - 1] != '/');

          UtilCreateCleanDirectory(a_directory_name, true);
        
          FAIOHeader<T> header(a_data);
          string filename = a_directory_name + string("/headerfile");
          std::ofstream headerfile;
          headerfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
          headerfile << header << std::endl;
          headerfile.flush();
          headerfile.close();

//        int coordinatorProc(ParallelDescriptor::IOProcessorNumber());
//        long bytesWritten(0);
//
//        FabArrayIO<T>::Header hdr(mf);
//
//        std::string filePrefix(mf_name + FabFileSuffix);
//
//        NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf);
//
//        for( ; nfi.ReadyToWrite(); ++nfi) 
//        {
//	  // ---- find the total number of bytes including fab headers if needed
//          const FABio &fio = FArrayBox::getFABio();
//          int whichRDBytes(whichRD->numBytes()), nFABs(0);
//          long writeDataItems(0), writeDataSize(0);
//          for(MFIter mfi(mf); mfi.isValid(); ++mfi) 
//          {
//	    const FArrayBox &fab = mf[mfi];
//	    bytesWritten += fab.box().numPts() * mf.nComp() * whichRDBytes;
//	    ++nFABs;
//	  }
//	  char *allFabData(nullptr);
//	  bool canCombineFABs(false);
//
//          allFabData = new(std::nothrow) char[bytesWritten];
//
//          long writePosition(0);
//          for(MFIter mfi(mf); mfi.isValid(); ++mfi) 
//          {
//            int hLength(0);
//            const FArrayBox &fab = mf[mfi];
//            writeDataItems = fab.box().numPts() * mf.nComp();
//            writeDataSize = writeDataItems * whichRDBytes;
//            char *afPtr = allFabData + writePosition;
//              
//            memcpy(afPtr + hLength, fab.dataPtr(), writeDataSize);
//
//            writePosition += hLength + writeDataSize;
//          }
//          nfi.Stream().write(allFabData, bytesWritten);
//          nfi.Stream().flush();
//          delete [] allFabData;
//        }
//
//
//        VisMF::FindOffsets(mf, filePrefix, hdr, groupSets, currentVersion,
//                           useDynamicSetSelection, nfi);
//
//        bytesWritten += VisMF::WriteHeader(mf_name, hdr, coordinatorProc);
//
//        delete whichRD;
//
//        return bytesWritten;
        }

    };
}
#endif
