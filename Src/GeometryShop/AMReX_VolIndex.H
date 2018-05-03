#ifndef _VOLINDEX_H_
#define _VOLINDEX_H_

#include "AMReX_REAL.H"
#include "AMReX_IntVect.H"
#include "AMReX_BaseIndex.H"

namespace amrex
{
  /// Volume of Fluid Index
  /**
     VolIndex is a very lightweight object
     to be used to distinguish between volumes of
     fluid in a cell.  We use VolIndex for
     array indexing.   Non-trivial constructors
     are all private because only EBISBox
     is supposed to make valid volume indicies.
  */
  class VolIndex: public BaseIndex
  {
  public:
    ///
    VolIndex();

    virtual ~VolIndex();

    ///
    //VolIndex& operator= (const VolIndex& a_vofin);

    ///
    
    bool operator== (const VolIndex& a_vofin) const {
	return ((m_cellIndex == a_vofin.m_cellIndex) &&
		(m_iv == a_vofin.m_iv));
    }

    ///
    
    const IntVect& gridIndex() const;

    ///
    
    int cellIndex() const;

    ///
    
    bool isDefined() const;

    ///
    // VolIndex(const VolIndex& a_vofin);

    ///
    void define(const VolIndex& a_vofin);

    ///
    
    VolIndex(const IntVect& a_ix,const  int& a_vofID);

    ///
    
    void define(const IntVect& a_ix,const  int& a_vofID);

    
    bool operator<(const VolIndex& rhs) const;

    
    /*@{*/
    int linearSize() const;
    void linearOut(void* const a_outBuf) const;
    void linearIn(const void* const inBuf);
    /*@}*/

    bool operator!=(const VolIndex& rhs) const;

    ///for stencil gymnastics
    void shift(const IntVect& a_iv)
    {
      m_iv += a_iv;
    }

  private:

    //cell of the VolIndex
    IntVect m_iv;

    //which VoF in the cell (0 to nVoF-1)
    int m_cellIndex;

    bool m_isDefined;
  };

  std::ostream& operator<< (std::ostream&       os,
                            const VolIndex& iv);

  inline bool VolIndex::operator<(const VolIndex& rhs) const
  {
    if (m_iv < (rhs.m_iv))
    {
      return true;
    }
    if (m_iv == rhs.m_iv)
    {
      return m_cellIndex < rhs.m_cellIndex;
    }
    return false;
  }

}
#endif
