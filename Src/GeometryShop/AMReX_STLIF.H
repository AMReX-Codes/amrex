
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

#ifndef _STLIF_H_
#define _STLIF_H_

#include <iostream>
#include <fstream>
using std::istream;
using std::ifstream;
using std::string;
using std::shared_ptr;

#include "AMReX_BaseIF.H"
#include "AMReX_STLExplorer.H"

namespace amrex
{
///
/**
   This implicit function reads an STL file and uses the polygonal information
   to provide edge intersections.  As such, calling it's "value" function is
   an error.  It is handled specially in "GeometryShop".  This is hardwired for ascii data.
*/
  class STLIF: public BaseIF
  {
  public:

    ///
    /**
       Constructor specifying filename (a_filename), the form of the data
       (a_dataType - ASCII or Binary), level set value (a_value), and whether
       inside the domain is where data is less than the level set value
       (a_inside).  Data is read from the file named and a complete ASCII
       header (see above) is expected.
    */
    STLIF(const string&     a_filename);


    /// Copy constructor
    STLIF(const STLIF& a_inputIF);

    /// Destructor
    virtual ~STLIF();

    ///
    /**
       For STLIF, calling this method is an error.
    */
    virtual Real value(const RealVect& a_point) const;

    virtual BaseIF* newImplicitFunction() const;

    shared_ptr<STLExplorer> getExplorer() const
      {
        return m_explorer;
      }

  protected:
    void makeExplorer();

    string                  m_filename;
    shared_ptr<STLExplorer> m_explorer;

  private:
    STLIF()
      {
      }

    void operator=(const STLIF& a_inputIF)
      {
      }
  };
}

#endif
