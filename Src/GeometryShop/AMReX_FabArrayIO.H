#ifndef AMREX_FabArrayIO_H_
#define AMREX_FabArrayIO_H_


#include <cmath>
#include <cstdlib>

#include "AMReX_VisMF.H"
#include "AMReX_EBArith.H"
#include "AMReX_SPMD.H"
#include "AMReX_Print.H"
#include "AMReX_Utility.H"
#include "AMReX_parstream.H"
#include "AMReX_BoxArray.H"


namespace amrex
{
  ///
  class FAIOElement
  {
  public:

    //every proc has its own file
    static string getFilename(const int& a_procid);

    ///
    FAIOElement (int a_fileid, long a_head, int a_boxid, long a_boxlen);
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
    void linearIn(const void* const buffer );

    string      m_filename;
    int         m_procid;  // used to create the filename
    long        m_head;    // Offset to start of FAB in file.
    int         m_boxid;   // integer into the box array
    long        m_boxlen;  //size of this box's data
  };


  ///
  class FAIOElementVec
  {
  public:


    ///
    FAIOElementVec(const Vector<FAIOElement>& a_vec)
      {
        m_vec = a_vec;
      }

    FAIOElementVec()
      {
      }

    ///
    int linearSize() const
      {
        return linearListSize(m_vec);
      }

    ///
    void linearOut(void* buffer ) const
      {
        linearListOut(buffer, m_vec);
      }

    ///
    void linearIn(const void* const buffer )
      {
        linearListIn(m_vec, buffer);
      }

    ///
    Vector<FAIOElement> m_vec;

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

    ///for writes, need to be able define the header
    FAIOHeader(const FabArray<T> & a_data)
      {
        m_ba = a_data.boxArray();
        Vector<FAIOElement> localElements(a_data.local_size());
        int ielem  = 0;
        int offset = 0;

        m_nComp = a_data.nComp();
        for(MFIter mfi(a_data); mfi.isValid(); ++mfi) 
        {
          int boxid   = mfi.index();
          //writes one file per proc
          //anything more flexible requires a bunch of communication
          int filenum = ParallelDescriptor::MyProc();
          long numBytesThisBox = a_data[mfi].nBytesFull();
          localElements[ielem] = FAIOElement(filenum, offset, boxid, numBytesThisBox);
          ielem++;
          offset += numBytesThisBox;
        }
        FAIOElementVec localVec(localElements);
        Vector<FAIOElementVec> allVecs;
        gather(allVecs, localVec, 0);
        FAIOElementVec globalVec;
        if(ParallelDescriptor::MyProc() == 0)
        {
          m_vecElem.resize(0);
          for(int ivec = 0; ivec < allVecs.size();  ivec++)
          {
            const Vector<FAIOElement>& onevec = allVecs[ivec].m_vec;
            globalVec.m_vec.insert(globalVec.m_vec.end(), onevec.begin(), onevec.end());
          }
        }
        broadcast(globalVec, 0);
        m_vecElem = globalVec.m_vec;
        std::sort(m_vecElem.begin(), m_vecElem.end()); //this should get them in the same order as the boxarray
      }

    ///For reads, need null construction
    FAIOHeader()
      {
      }

    BoxArray                   m_ba;        
    Vector<FAIOElement>   m_vecElem;   
    int m_nComp;
  private:

  };

  ///write a header to   disk in ascii
  template <class T>
  std::ostream& 
  operator<< (std::ostream& a_os, FAIOHeader<T>& a_hd)
  {
    a_os << a_hd.m_nComp      << endl;
    a_os << a_hd.m_vecElem.size() << endl;
    for(int ivec = 0; ivec < a_hd.m_vecElem.size(); ivec++)
    {
      a_os << a_hd.m_vecElem[ivec] << endl;
    }
    long bxsize = a_hd.m_ba.size();

    //box array does not have proper ostream, istream functions.
    a_os << bxsize << endl;
    for(int ibox = 0; ibox < bxsize; ibox++)
    {
      a_os << a_hd.m_ba[ibox] << endl;
    }
    return a_os;
  }

  ///read a header from disk in ascii
  template <class T>
  std::istream& 
  operator>> (std::istream& a_is, FAIOHeader<T>& a_hd)
  {
    a_is >> a_hd.m_nComp;
    int isize;
    a_is >> isize;
    a_hd.m_vecElem.resize(isize);
    for(int ivec = 0; ivec < a_hd.m_vecElem.size(); ivec++)
    {
      a_is >> a_hd.m_vecElem[ivec];
    }

    long bxsize;

    //box array does not have proper ostream, istream functions.
    a_is >> bxsize;
    Vector<Box> boxes(bxsize);
    for(int ibox = 0; ibox < bxsize; ibox++)
    {
      a_is >> boxes[ibox];
    }
    
    a_hd.m_ba = BoxArray(&boxes[0], bxsize);
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
    /// write to disk.  one file per proc
    static void
    write(const FabArray<T>&    a_data,
          const std::string&    a_directory_name)
      {
        BL_PROFILE("FabArrayIO::write");
        BL_ASSERT(a_directory_name[a_directory_name.length() - 1] != '/');

        UtilCreateCleanDirectory(a_directory_name, true);

        //first make and write the header.  This is the only blocking bit of this
        FAIOHeader<T> header(a_data);
        string filename = a_directory_name + string("/headerfile");
        std::ofstream headerfile;
        headerfile.open(filename.c_str(), std::ios::out | std::ios::trunc);
        headerfile << header << std::endl;
        headerfile.flush();
        headerfile.close();

        //now allocate all the big string of bytes to write to disk. 
        //first we need to find the size.  Since every proc writes to its 
        //own file, we only need to allocate once and write once and there
        // is no communication here
        size_t totallength = 0;

        for(MFIter mfi(a_data); mfi.isValid(); ++mfi) 
        {
          size_t boxlen = a_data[mfi].nBytesFull();
          totallength += boxlen;
        }
        char* outbuf = new char[totallength];
        string datafilename = a_directory_name + string("/") + FAIOElement::getFilename(ParallelDescriptor::MyProc());

        //now do the actual write.  Again, since every proc as its own file,
        //I can do a simple open and close of the file
        std::fstream fileStream;
        fileStream.open(datafilename.c_str(),
                        std::ios::out | std::ios::trunc | std::ios::binary);
          
        //no need for fancy seeks. this ought to work just
        //putting one in after the next
        char* movingBuf = outbuf;
        for(MFIter mfi(a_data); mfi.isValid(); ++mfi) 
        {
          size_t boxlen = a_data[mfi].nBytesFull();
          a_data[mfi].copyToMemFull(movingBuf);
          //pout() << "write: filename = " << datafilename << ", box = " << mfi.validbox() << endl;
          movingBuf += boxlen;
        }
        fileStream.write(outbuf, totallength);
        fileStream.flush();
        fileStream.close();

        delete[] outbuf;
      }


    ///read from disk --null dm means make my own
    static void
    read(FabArray<T>      &    a_data,
         const std::string&    a_directory_name,
         const DistributionMapping * dmPtr = NULL)
      {
        BL_PROFILE("FabArrayIO::read");
        //first make and write the header.  This is the only blocking bit of this
        FAIOHeader<T> header;
        string filename = a_directory_name + string("/headerfile");
        std::ifstream headerfile;
        headerfile.open(filename.c_str(), std::ios::in);
        headerfile >> header;
        headerfile.close();

        //define with default factory and have the I/O fill in the holders
        DistributionMapping dm;
        if(dmPtr == NULL)
        {
          dm = DistributionMapping(header.m_ba);
        }
        else
        {
          dm = *dmPtr;
        }

        a_data.define(header.m_ba, dm,  header.m_nComp, 0, MFInfo(), DefaultFabFactory<T>());
        int ibox = 0;
        for(MFIter mfi(a_data); mfi.isValid(); ++mfi) 
        {
          int baindex = mfi.index();
          const FAIOElement& elem = header.m_vecElem[baindex];
          string datafile = a_directory_name + "/" + elem.m_filename;
          //cannot do this one in one fell swoop because distribution may have changed
          std::fstream fileStream;
          fileStream.open(datafile.c_str(),          
                          std::ios::in  | std::ios::binary);

          //pout() << "read: filename = " << datafile << ", box = " << mfi.validbox() << endl;

          fileStream.seekg(elem.m_head, std::ios::beg);
          size_t boxsize = elem.m_boxlen;
          char*  inbuf = new char[boxsize];
          fileStream.read(inbuf, boxsize);

          a_data[mfi].copyFromMemFull(inbuf);
          fileStream.flush();
          fileStream.close();
          delete[] inbuf;
          ibox++;
        }
                      
      }

  };
}
#endif
