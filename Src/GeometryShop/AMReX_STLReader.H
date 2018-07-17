#ifndef AMREX_STLREADER_H_
#define AMREX_STLREADER_H_

#include <iostream>
#include <fstream>
using std::istream;
using std::ifstream;
#include "AMReX_STLMesh.H"

///base class so that the reader can have the generality of reading binary or ASCII STL data.
namespace amrex
{
  class STLReader
  {
  public:
    /// Constructor - read from standard input
    STLReader() {}

    /// Destructor
    virtual ~STLReader() {}

    /// Return header information
    virtual string* GetHeader() const = 0;

    /// Return number of elements
    virtual void GetNtri(int& a_ntri) const = 0;

    /// Return whether number of elements from header matches file
    virtual void GetNtriMatch(bool& a_ntriMatch) const = 0;

    /// Return pointer to the mesh
    virtual shared_ptr<STLMesh> GetMesh() const = 0;

  protected:
    virtual void ReadData(istream& a_file,
                          const int offset) = 0;

    int  m_ntri;      // number of triangles read in
    bool m_ntriMatch; // true if m_ntri equals the number of triangles declared in the header

    // actual data - shared by all copies
    shared_ptr<STLMesh> m_stlmesh; // pointer to the mesh

  private:
    void operator=(const STLReader& a_inputReader)
      {
      }
  };
}

#endif
