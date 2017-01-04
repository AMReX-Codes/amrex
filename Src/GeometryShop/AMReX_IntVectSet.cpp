

/*
 *      .o.       ooo        ooooo ooooooooo.             ooooooo  ooooo 
 *     .888.      `88.       .888' `888   `Y88.            `8888    d8'  
 *    .8"888.      888b     d'888   888   .d88'  .ooooo.     Y888..8P    
 *   .8' `888.     8 Y88. .P  888   888ooo88P'  d88' `88b     `8888'     
 *  .88ooo8888.    8  `888'   888   888`88b.    888ooo888    .8PY888.    
 * .8'     `888.   8    Y     888   888  `88b.  888    .o   d8'  `888b   
 *o88o     o8888o o8o        o888o o888o  o888o `Y8bod8P' o888o  o88888o 
 *
 */

#include "AMReX_Box.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_IntVectSet.H"
namespace amrex
{
  ///
  IntVectSet::
  IntVectSet(const Box& a_box)
  {
    for(BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        *this |= bit();
      }
  }
  ///
  IntVectSet::
  IntVectSet(const IntVectSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
  ///
  void 
  IntVectSet::
  define(const Box& a_box)
  {
    *this = IntVectSet(a_box);
  }
  ///
  void 
  IntVectSet::
  define(const IntVectSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
  }
  ///
  IntVectSet& 
  IntVectSet::
  operator=(const IntVectSet& a_sivs)
  {
    m_stdSet = a_sivs.m_stdSet; 
    return *this;
  }
  ///
  IntVectSet& 
  IntVectSet::
  operator|=(const IntVectSet& a_sivs)
  {
    const std::set<IntVect, lex_compare_iv> inputset = a_sivs.m_stdSet;
    std::set<IntVect,lex_compare_iv>::iterator it;
    for(it = inputset.begin(); it!=  inputset.end(); ++it)
      {
        m_stdSet.insert(*it);
      }
    return *this;
  }
  ///
  IntVectSet & 
  IntVectSet::
  operator|=(const IntVect& a_iv)
  {
    m_stdSet.insert(a_iv);
    return *this;
  }
  ///
  IntVectSet& 
  IntVectSet::
  operator|=(const Box& a_box)
  {
    for(BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        m_stdSet.insert(bit());
      }
    return *this;
  }
  ///
  IntVectSet& 
  IntVectSet::
  operator&=(const IntVectSet& a_sivs)
  {
    if(&a_sivs != this)
      {
        std::set<IntVect, lex_compare_iv> newSet;
        std::set<IntVect, lex_compare_iv>::iterator it;
        for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
          {
            const IntVect& iv = *it;
            if(contains(iv) && a_sivs.contains(iv))
              {
                newSet.insert(iv);
              }
          }
        m_stdSet = newSet;
      }
    return *this;
  }

  ///and
  IntVectSet& 
  IntVectSet::
  operator&=(const Box& a_box)
  {
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        if(!a_box.contains(iv))
          {
            m_stdSet.erase(it);
          }
      }
    return *this;
  }
  ///not
  IntVectSet& 
  IntVectSet::
  operator-=(const IntVectSet& a_sivs)
  {
    std::set<IntVect, lex_compare_iv>::iterator it;
    //leaving out the ++it because  erase 
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); )
      {
        if(a_sivs.contains(*it))
          {
            m_stdSet.erase(it++);
          }
        else
          {
            ++it;
          }
      }
    return *this;
  }
  ///not
  IntVectSet& 
  IntVectSet::
  operator-=(const IntVect& a_iv)
  {
    if(contains(a_iv))
      {
        m_stdSet.erase(m_stdSet.find(a_iv));
      }
    return *this;
  }
  ///not
  IntVectSet& 
  IntVectSet::
  operator-=(const Box& a_box)
  {
    for(BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        *this -= bit();
      }
    return *this;
  }
  ///
  bool 
  IntVectSet::
  operator==(const IntVectSet& a_lhs) const
  {
    if(a_lhs.m_stdSet.size() != m_stdSet.size())
      {
        return false;
      }

    bool retval = true;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it!=  m_stdSet.end(); ++it)
      {
        if((!contains(*it)) || (!a_lhs.contains(*it)))
          {
            retval = false;
            break;
          }
      }
    return retval;
  }

  ///
  bool 
  IntVectSet::
  contains(const IntVect& a_iv) const
  {
    std::set<IntVect, lex_compare_iv>::iterator it = m_stdSet.find(a_iv);
    return (it != m_stdSet.end());
  }

  ///
  bool 
  IntVectSet::
  contains(const Box& a_box) const
  {
    bool retval = true;
    for(BoxIterator bit(a_box); bit.ok(); ++bit)
      {
        if(!contains(bit()))
          {
            retval = false;
            break;
          }
      }
    return retval;
  }

  ///
  void 
  IntVectSet::
  grow(int igrow)
  {
    IntVectSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        Box grid(iv, iv);
        grid.grow(igrow);
        newSet |= grid;
      }
    *this = newSet;
  }

  ///
  void 
  IntVectSet::
  grow(int idir, int igrow)
  {
    IntVectSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        Box grid(iv, iv);
        grid.grow(idir, igrow);
        newSet |= grid;
      }
    *this = newSet;
  }

  ///
  void 
  IntVectSet::
  growHi()
  {
    for(int idir = 0; idir < SpaceDim; idir++)
      {
        growHi(idir);
      }
  }

  ///
  void 
  IntVectSet::
  growHi(int a_dir)
  {
    IntVectSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        Box grid(iv, iv);
        grid.growHi(a_dir);
        newSet |= grid;
      }
    *this = newSet;
  }

  ///
  void 
  IntVectSet::
  refine(int iref)
  {
    IntVectSet newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        Box grid(iv, iv);
        grid.refine(iref);
        newSet |= grid;
      }
    *this = newSet;
  }

  ///
  void 
  IntVectSet::
  coarsen(int iref)
  {
    std::set<IntVect, lex_compare_iv> newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        Box grid(iv, iv);
        grid.coarsen(iref);
        for(BoxIterator bit(grid); bit.ok(); ++bit)
          {
            newSet.insert(bit());
          }
      }
    m_stdSet = newSet;
  }

  ///
  void 
  IntVectSet::
  shift(const IntVect& a_iv)
  {
    std::set<IntVect, lex_compare_iv> newSet;
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        IntVect iv = *it;
        iv.shift(a_iv);
        newSet.insert(iv);
      }
    m_stdSet = newSet;
  }

  ///
  void 
  IntVectSet::
  clear()
  {
    std::set<IntVect, lex_compare_iv> newSet;
    m_stdSet = newSet;
  }

  ///
  Box 
  IntVectSet::
  minBox() const
  {
    int bignum = 100000;
    IntVect lo = bignum*IntVect::TheUnitVector();
    IntVect hi =-bignum*IntVect::TheUnitVector();
    std::set<IntVect, lex_compare_iv>::iterator it;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        const IntVect& iv = *it;
        for(int idir = 0; idir < SpaceDim; idir++)
          {
            lo[idir] = std::min(lo[idir], iv[idir]);
            hi[idir] = std::max(hi[idir], iv[idir]);
          }
      }
  
    Box retval(lo, hi);
    return retval;
  }

  ///
  bool 
  IntVectSet::
  isEmpty() const
  {
    return (m_stdSet.size() == 0);
  }
  ///
  void 
  IntVectSet::
  getVectorIV(std::vector<IntVect>& a_vect) const
  {
    a_vect.resize(m_stdSet.size());

    std::set<IntVect, lex_compare_iv>::iterator it;
    int ivec = 0;
    for(it = m_stdSet.begin(); it != m_stdSet.end(); ++it)
      {
        a_vect[ivec] = *it;
        ivec++;
      }
  }
  ///
  void 
  IntVectSet::
  makeEmpty() 
  {
    clear();
  }

  ///
  int
  IntVectSet::
  numPts() const
  {
    return m_stdSet.size();
  }

  ///
  void 
  IntVectSet::
  define(const std::vector<IntVect>& a_vect)
  {
    makeEmpty();
    for(int ivec = 0; ivec  < a_vect.size(); ivec++)
      {
        (*this ) |= a_vect[ivec];
      }
  }



  IVSIterator::
  IVSIterator()
  {
    m_ivs = NULL;
  }

  ///
  IVSIterator::
  IVSIterator(const IntVectSet& ivs)
  {
    m_ivs = &ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

  ///
  void 
  IVSIterator::
  define(const IntVectSet & a_ivs)
  {
    m_ivs = &a_ivs;
    m_iter = m_ivs->m_stdSet.begin();
  }

  ///
  const IntVect& 
  IVSIterator::
  operator()() const 
  {
    return *m_iter;
  }

  ///
  bool 
  IVSIterator::
  ok() const
  {
    assert(m_ivs != NULL);
    return (m_iter != m_ivs->m_stdSet.end());
  }

  ///
  void 
  IVSIterator::
  operator++()
  {
    m_iter++;
  }

  ///
  void 
  IVSIterator::
  begin()
  {
    assert(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.begin();
  }

  ///
  void 
  IVSIterator::
  end()
  {
    assert(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.end();
  }

  ///
  void 
  IVSIterator::
  clear()
  {
    m_ivs = NULL;
  }
}
