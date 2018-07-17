#include <cmath>
#include <cstdlib>

#include "AMReX_FabArrayIO.H"

namespace amrex
{
  string
  FAIOElement::
  getFilename(const int& a_fileid)
  {
    string retval = string("data_") + EBArith::convertInt(a_fileid);
    return retval;
  }

  ///
  FAIOElement::
  FAIOElement (int a_fileid, long a_head, int a_boxid, long a_boxlen)
  {
    m_procid = a_fileid;
    m_head   = a_head;
    m_boxid  = a_boxid;
    m_boxlen = a_boxlen;
    m_filename = getFilename(a_fileid);
  }

  ///
  bool
  FAIOElement::
  operator<(const FAIOElement& a_fab) const
  {
    return m_boxid < a_fab.m_boxid;
  }

  ///
  int 
  FAIOElement::
  linearSize() const
  {
    int retval = 0;
    //procid
    retval += sizeof(int);
    //boxid
    retval += sizeof(int);
    //head
    retval += sizeof(long);
    //boxlen
    retval += sizeof(long);
    return retval;
  }

  ///
  void 
  FAIOElement::
  linearOut(void* buffer ) const
  {
    int* intbuf = (int*)buffer;
    *intbuf = m_procid; 
    intbuf++;
    *intbuf = m_boxid;
    intbuf++;
    long* longbuf = (long*)(intbuf);
    *longbuf = m_head;
    longbuf++;
    *longbuf = m_boxlen;
    longbuf++;
  }

  ///
  void 
  FAIOElement::
  linearIn(const void* const buffer )
  {
    int* intbuf = (int*)buffer;
    m_procid = *intbuf;
    intbuf++;
    m_boxid = *intbuf;
    intbuf++;
    long* longbuf = (long*)(intbuf);
    m_head = *longbuf;
    longbuf++;
    m_boxlen = *longbuf;
    longbuf++;
    m_filename = getFilename(m_procid);
  }

  ///write a header to   disk in ascii
  std::ostream& 
  operator<< (std::ostream& a_os, FAIOElement& a_elem)
  {
    a_os << a_elem.m_filename << "   ";
    a_os << a_elem.m_procid   << "   ";
    a_os << a_elem.m_head     << "   ";
    a_os << a_elem.m_boxid    << "   ";
    a_os << a_elem.m_boxlen    << "   ";
    return a_os;
  }

  ///read a header from disk in ascii
  std::istream& 
  operator>> (std::istream& a_is, FAIOElement& a_elem)
  {
    a_is >> a_elem.m_filename;
    a_is >> a_elem.m_procid  ;
    a_is >> a_elem.m_head    ;
    a_is >> a_elem.m_boxid   ;
    a_is >> a_elem.m_boxlen   ;
    return a_is;
  }
}

