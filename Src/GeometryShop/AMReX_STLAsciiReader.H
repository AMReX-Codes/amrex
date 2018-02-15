
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

#ifndef _STLASCIIREADER_H_
#define _STLASCIIREADER_H_

#include <iostream>
#include <string>
using namespace std;


#include "STLReader.H"
#include "STLMesh.H"


namespace amrex
{
///
/**
 * Reads ascii STL files and generates a mesh
 * see http://www.cplusplus.com/doc/tutorial/files/
 */
  class STLAsciiReader: public STLReader
  {
  public:
    /// Constructor - read from standard input
    STLAsciiReader();

    /// Constructor - read from file name
    STLAsciiReader(const string& a_filename);

    /// Destructor
    ~STLAsciiReader();

    /// Return header information
    string* GetHeader() const;

    /// Return number of elements
    void GetNtri(int& a_ntri) const;

    /// Return whether number of elements from header matches file
    void GetNtriMatch(bool& a_ntriMatch) const;

    /// Return pointer to the mesh
    RefCountedPtr<STLMesh> GetMesh() const;

  protected:
    void ReadData(istream&   a_file,
                  const int offset);

    string* m_header;    // header info
    int     m_ntri;      // number of triangles read in
    bool    m_ntriMatch; // true if m_ntri equals the number of triangles declared in the header

    // actual data - shared by all copies
    RefCountedPtr<STLMesh> m_stlmesh; // pointer to the mesh

  private:
    void operator=(const STLAsciiReader& a_inputReader)
      {
      }
  };

}
#endif
