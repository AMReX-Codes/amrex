#include "AMReX_Box.H"
#include "AMReX_BoxIterator.H"
#include "AMReX_IntVectSet.H"

namespace amrex
{
  ///
  size_t 
  IntVectSet::
  linearSize() const
  {
    //the one is for the size
    size_t retval = sizeof(int);
    retval += numPts()*IntVect::linearSize();
    return retval;
  }

  ///
  void 
  IntVectSet::
  linearOut(void* a_buffer ) const
  {
    int* intbuf = (int *) a_buffer;
    int numpts = numPts();
    *intbuf = numpts;
    intbuf++;
    unsigned char* buf = (unsigned char*) intbuf;
    int icount = 0;
    for(IVSIterator ivsit(*this); ivsit.ok(); ++ivsit)
    {
      const IntVect& iv = ivsit();
      iv.linearOut(buf);
      buf += IntVect::linearSize();
      icount++;
    }
    BL_ASSERT(icount == numpts);
  }

  ///
  void 
  IntVectSet::
  linearIn(const void* const a_buffer )
  {
    int* intbuf = (int *) a_buffer;
    int numpts = *intbuf;
    intbuf++;
    unsigned char* buf = (unsigned char*) intbuf;
    for(int ivec = 0; ivec < numpts; ivec++)
    {
      IntVect iv;
      iv.linearIn(buf);
      buf += IntVect::linearSize();

      *this |= iv;
    }
  }

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
    m_stdSet.insert(a_sivs.cbegin(), a_sivs.cend());
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
        std::unordered_set<IntVect, IntVect::shift_hasher> newSet;
	for (const auto& iv : *this) {
	    if (a_sivs.contains(iv)) {
		newSet.insert(iv);
	    }
	}
	m_stdSet = std::move(newSet);
      }
    return *this;
  }

  ///and
  IntVectSet& 
  IntVectSet::
  operator&=(const Box& a_box)
  {
      *this = intersect(*this,a_box);
      return *this;
  }

  ///not
  IntVectSet& 
  IntVectSet::
  operator-=(const IntVectSet& a_sivs)
  {
    //leaving out the ++it because  erase 
    for(auto it = m_stdSet.begin(); it!=  m_stdSet.end(); )
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
    m_stdSet.erase(a_iv);
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
      return this->m_stdSet == a_lhs.m_stdSet;
  }

  ///
  bool 
  IntVectSet::
  contains(const IntVect& a_iv) const
  {
    auto it = m_stdSet.find(a_iv);
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
  bool 
  IntVectSet::
  contains(const IntVectSet& a_ivs) const
  {
    bool retval = true;
    for(IVSIterator ivsit(a_ivs); ivsit.ok(); ++ivsit)
      {
        if(!contains(ivsit()))
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
    for(const auto& iv : *this)
      {
        Box grid(iv, iv);
        grid.grow(igrow);
        newSet |= grid;
      }
    *this = std::move(newSet);
  }

  ///
  void 
  IntVectSet::
  grow(int idir, int igrow)
  {
    IntVectSet newSet;
    for(const auto& iv : *this)
      {
        Box grid(iv, iv);
        grid.grow(idir, igrow);
        newSet |= grid;
      }
    *this = std::move(newSet);
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
    for(const auto& iv : *this)
      {
        Box grid(iv, iv);
        grid.growHi(a_dir);
        newSet |= grid;
      }
    *this = std::move(newSet);
  }

  ///
  void 
  IntVectSet::
  refine(int iref)
  {
    IntVectSet newSet;
    for(const auto& iv : *this)
      {
        Box grid(iv, iv);
        grid.refine(iref);
        newSet |= grid;
      }
    *this = std::move(newSet);
  }

  ///
  void 
  IntVectSet::
  coarsen(int iref)
  {
    container_type newSet;
    for(auto ivcoar : *this)
      {
        ivcoar.coarsen(iref);
        newSet.insert(ivcoar);
      }
    m_stdSet = std::move(newSet);
  }

  ///
  void 
  IntVectSet::
  shift(const IntVect& a_iv)
  {
    container_type newSet;
    for(auto iv : *this)
      {
        iv.shift(a_iv);
        newSet.insert(iv);
      }
    m_stdSet = std::move(newSet);
  }

  ///
  void 
  IntVectSet::
  clear()
  {
    m_stdSet.clear();
  }

  ///
  Box 
  IntVectSet::
  minBox() const
  {
    int bignum = 100000;
    IntVect lo = bignum*IntVect::TheUnitVector();
    IntVect hi =-bignum*IntVect::TheUnitVector();
    for (const auto& iv : *this)
      {
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
      return m_stdSet.empty();
  }
  ///
  void 
  IntVectSet::
  getVectorIV(Vector<IntVect>& a_vect) const
  {
    a_vect.resize(m_stdSet.size());
    std::copy(std::begin(*this), std::end(*this), std::begin(a_vect));
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
  define(const Vector<IntVect>& a_vect)
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
    BL_ASSERT(m_ivs != NULL);
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
    BL_ASSERT(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.begin();
  }

  ///
  void 
  IVSIterator::
  end()
  {
    BL_ASSERT(m_ivs != NULL);
    m_iter = m_ivs->m_stdSet.end();
  }

  ///
  void 
  IVSIterator::
  clear()
  {
    m_ivs = NULL;
  }

    IntVectSet intersect (const IntVectSet& a, const IntVectSet& b, const Box& bx)
    {
	IntVectSet r;
	for (const auto& iv : a) {
	    if (bx.contains(iv) && b.contains(iv)) {
		r |= iv;
	    }
	}
	return r;
    }

    IntVectSet intersect (const IntVectSet& ivs, const Box& bx)
    {
	IntVectSet r;
	for (const auto& iv : ivs) {
	    if (bx.contains(iv)) {
		r |= iv;
	    }
	}
	return r;
    }
}
