
#ifndef AMREX_FNERGINDEXEDMOMENTS_H_
#define AMREX_FNERGINDEXEDMOMENTS_H_

#include "AMReX_IndexTM.H"
#include "AMReX_Vector.H"
#include "AMReX_SPACE.H"
#include "AMReX_MonomialPowers.H"
#include "AMReX_MomentIterator.H"
#include "AMReX_Factorial.H"
#include <map>

using std::map;
using amrex::Vector;

//! \class IndexedMoments
//! Vector-like container for multi-indexed Real values up to 
//! some max multi-index P (sum of indicies <= P). 
//! Layout is 0th dim first, then 1st, etc.
//! \tparam D The Dimension of the container
template <int Dim, int P>
class IndexedMoments
{
public:
  /// constructor---make statics first time called
  IndexedMoments()
    {
      m_moms.resize(size());
      m_isRegular = 0;
    }

  //! Destructor
  ~IndexedMoments() {};

  static int size()
    {
      return MonomialPowers<Dim,P>::size();
    }

  //! Retrieve the moment from the index 
  //! \params a_index The lexigraphical index to the data
  Real& operator[](int a_i) 
    {
      return m_moms[a_i]; 
    };

  //! Retrieve the moment from the index 
  //! \params a_index The lexigraphical index to the data
  const Real& operator[](int a_i)  const
    {
      return m_moms[a_i]; 
    };

  //! Retrieve the moment from the index 
  //! \params a_index The multi-index that's needed
  Real& operator[](const IndexTM<int,Dim>& a_index) 
    {
      return m_moms[indexOf(a_index)]; 
    };

  //! Retrieve the moment from the index 
  //! \params a_index The multi-index that's needed
  const Real& operator[](const IndexTM<int,Dim>& a_index) const
    {
      return m_moms[indexOf(a_index)]; 
    };

  ///add compoenentwise
  IndexedMoments<Dim, P>&  operator+=(const IndexedMoments<Dim, P>& a_input)
    {
      for(MomentIterator<Dim,P> it; it.ok(); ++it)
      {
        (*this)[it()] += a_input[it()];
      }
      return *this;
    }

  ///multiply each component by constant
  IndexedMoments<Dim, P>&  operator*=(const Real& a_factor)
    {
      for(MomentIterator<Dim,P> it; it.ok(); ++it)
      {
        (*this)[it()] *= a_factor;
      }
      return *this;
    }



  ///for linearization
  static size_t linearSize()
    {
      int numReals = size();
      //add in the real numbers
      size_t retval = numReals*sizeof(Real);

      //add in the m_isRegular
      retval += sizeof(int);
    
      return retval;
    }

  void linearOut(void* const a_outbuf) const
    {
      Real* realbuf = (Real*) a_outbuf;
      for(int ivec = 0; ivec < size(); ivec++)
      {
        *realbuf = m_moms[ivec];
        ++realbuf;
      }
      int* intbuf = (int *) realbuf;
      *intbuf =  m_isRegular;
    }

  void linearIn(const void* const a_inbuf)
    {
      Real* realbuf = (Real*) a_inbuf;
      for(int ivec = 0; ivec < size(); ivec++)
      {
        m_moms[ivec]= *realbuf;
        ++realbuf;
      }
      int* intbuf = (int *) realbuf;
      m_isRegular = *intbuf ;
    }

  /// set to a regular IndexTM
  void setRegular(const Real a_dx)
    {
      BL_ASSERT(a_dx > 0);
      m_isRegular = 1;

      // Get the set of multi-indexes in the correct order
      const Vector<IndexTM<int, Dim> >& mix = MonomialPowers<Dim,P>::getMonomialPowers();

      // pout() << "Regular moments:" << endl;
      for (int ix=0; ix < mix.size(); ++ix)
      {
        IndexTM<int,Dim> index = mix[ix];

        bool even = true;
        for (int d=0; (even) && (d < Dim); ++d)
          even = even && !(index[d] % 2);

        //! Only moments with an even multi-index are non-zero
        Real moment = 0;
        if (even)
        {
          moment = POW(0.5 * a_dx, index.sum()) * std::pow(a_dx, Dim);
          for (int d=0; d < Dim; ++d)
            moment /= (Real) (index[d]+1);
        }
        m_moms[ix] = moment;

        // pout() << ix << ": " << index << " = " << moment << endl;
      }

    }



  ///
  bool isRegular() const
    {
      return (m_isRegular==1);
    }

  /// for use with irregnode
  IndexedMoments<Dim, P>& operator=(const map<IndexTM<int, Dim>,  Real>& a_mapin)
    {
  
      for(MomentIterator<Dim, P> iter; iter.ok(); ++iter)
      {
        IndexTM<int, Dim> mono = iter();
        Real moment= getMoment(iter(), a_mapin);
        (*this)[iter()] = moment;
      }
      return *this;
    }


  /// 
  IndexedMoments<Dim, P>& operator=(const IndexedMoments<Dim, P>& a_input)
    {
      if(&a_input != this)
      {
        m_isRegular = a_input.m_isRegular;
        m_moms      = a_input.m_moms;
      }
      return *this;
    }
  

  ///
  /**
     shift moment by the input vector distance.
     this changes the current object from integral(x^p)
     to integral((x+x0)^p), where x0 = a_distance
  */
  void shift(const IndexTM<Real, Dim>& a_distance)
    {
      IndexTM<Real, Dim> x = a_distance;

      IndexedMoments<Dim, P> retval = *this;
      for(MomentIterator<Dim,P> itouter; itouter.ok(); ++itouter)
      {
        IndexTM<int, Dim> r = itouter();
        // - use the binomial shift rule to move them to a_x0
        Real moment = 0.0;
        for(MomentIterator<Dim,P> itinner; itinner.ok(); ++itinner)
        {
          IndexTM<int, Dim> q = itinner();
         
          if(q.componentwiseLE(r))
          {
            IndexTM<int, Dim> p = r - q;
            Real m0 = (*this)[q];
            Real xpow = 1;
            for(int idir = 0; idir < Dim; idir++)
            {
              xpow *=  std::pow(x[idir], p[idir]);
            }
              
            moment += pCk<Dim>(r, q) * xpow * m0;
            
          } //end if (q <= r)
        }//end loop over inner moments
        retval[r] = moment;
      }
      *this = retval;
    }
  
  ///
  void setToZero()
    {
      for(int ivec = 0; ivec < size(); ivec++)
      {
        m_moms[ivec] = 0.;
      }
    }

  /// Calculate what linear index this multi-index is
  static int indexOf(IndexTM<int,Dim> a_index)
    {
      BL_ASSERT(a_index.sum() <= P);
      BL_ASSERT(allPositive(a_index));

      // Formula for the layout index with 0th dim first, then 1st, etc.
      int index = a_index[0];
      int rem = P;
      for (int d=Dim-1; d > 0; --d)
      {
        index += s_mpowers_max_sizes[d][rem] - s_mpowers_max_sizes[d][rem-a_index[d]];
        rem -= a_index[d];
      }

      return index;
    }

  ///
  static IndexTM<int,Dim> getIndex(const int& a_linearIndex)
    {
      const Vector<IndexTM<int, Dim> >& multiIndices = MonomialPowers<Dim,P>::getMonomialPowers();
      return multiIndices[a_linearIndex];
    }

  ///
  /**
     Say <A> = sum_p(CA m^p),
     and <B> = sum_q(CB m^q).

     This sets the current data to 
     the set of coefficents M such that
     <AB> = sum_r(M m^r) + O(h^P+1).

     We drop all coefficents for powers s.t. p + q > P.
  */
  void setToTruncatedMultiply(const IndexedMoments<Dim, P> & a_CA,
                              const IndexedMoments<Dim, P> & a_CB)
    {
      this->setToZero();
      //here we are computing sum(m_p) = integral(ni(x-xbar)^p1/p! (x-xbar)^p2)
      // = sum_p2 ((n^p2)/p2!) m_p+p12
      for(MomentIterator<Dim, P> momitOuter; momitOuter.ok(); ++momitOuter)
      {
        //this is where the answer goes 
        IndexTM<int,BL_SPACEDIM>  po = momitOuter();

        for(MomentIterator<Dim, P> momitInner; momitInner.ok(); ++momitInner)
        {
          IndexTM<int,Dim> pi = momitInner();
          IndexTM<int,Dim> psum = po + pi;
          if(psum.sum() <= P)
          {
            Real incr = a_CA[pi]*a_CB[po];
            (*this)[psum] += incr;
          }
        }
      }

    }

  ///divides each entry by p!
  void divideByFactorial()
    {
      for(MomentIterator<Dim, P> momit; momit.ok(); ++momit)
      {
        IndexTM<int,Dim>  p = momit();
        Real factorial = pfactorial(p);
        (*this)[p] /= factorial;
      }
    }


  ///multiply each entry by p!
  void multiplyByFactorial()
    {
      for(MomentIterator<Dim, P> momit; momit.ok(); ++momit)
      {
        IndexTM<int,Dim>  p = momit();
        Real factorial = pfactorial(p);
        (*this)[p] *= factorial;
      }
}

protected:

  ///
  Real 
  getMoment(const IndexTM<int, Dim>        & a_mono,
            const map<IndexTM<int, Dim>,  Real>& a_mapin) const
    {
      /**/
      typename std::map<IndexTM<int, Dim>,  Real>::const_iterator iterinator;
      iterinator = a_mapin.find(a_mono);

      Real moment = 1.23456789e10;

      if (iterinator != a_mapin.end())
      {
        moment = iterinator->second;

      }
      else
      {
        std::cout << "monomial = " ;
        for(int idir = 0; idir < Dim; idir++)
        {
          std::cout << a_mono[idir] << " ";
        }
        std::cout << std::endl;
        amrex::Abort("getMoments in IndexedMoments: moment not found in map");
      }
      return moment;
    }

private:

  // Indicator that we contain only "full" moments
  int m_isRegular; 

  // Our moments to store    
  Vector<Real> m_moms;

}; 

/// Calculate what linear index this multi-index is 
///without the order stuff
template<int Dim>
int getIndMomLinearIndex(const IndexTM<int,Dim>& a_index, 
                         const int             & a_order)
{
  int retval= 0;
  if(a_order      == 1)
  {
    retval = IndexedMoments<Dim, 1>::indexOf(a_index);
  }                               
  else if(a_order == 2)             
  {                               
    retval = IndexedMoments<Dim, 2>::indexOf(a_index);
  }                               
  else if(a_order == 3)             
  {                               
    retval = IndexedMoments<Dim, 3>::indexOf(a_index);
  }                               
  else if(a_order == 4)             
  {                               
    retval = IndexedMoments<Dim, 4>::indexOf(a_index);
  }
  else
  {
    amrex::Error("need to add more cases to getIndMomLinearIndex");
  }
  return retval;
}

///
template<int Dim>
const IndexTM<int,Dim>
getIndMomMultiIndex(const int             & a_index,
                    const int             & a_order)
{
  IndexTM<int,Dim> retval;
  if(a_order      == 1)
  {
    retval = IndexedMoments<Dim,1>::getIndex(a_index);
  }
  else if(a_order == 2)
  {
    retval = IndexedMoments<Dim,2>::getIndex(a_index);
  }
  else if(a_order == 3)
  {
    retval = IndexedMoments<Dim,3>::getIndex(a_index);
  }
  else if(a_order == 4)
  {
    retval = IndexedMoments<Dim,4>::getIndex(a_index);
  }
  else
  {
    amrex::Error("need to add more cases to getIndMomMultiIndex");
  }
  return retval;
}

///
template<int Dim>
int
getIndMomSize(const int             & a_order)
{
  int retval = 0;
  if(a_order      == 1)
  {
    retval = IndexedMoments<Dim,1>::size();
  }
  else if(a_order == 2)
  {
    retval = IndexedMoments<Dim,2>::size();
  }
  else if(a_order == 3)
  {
    retval = IndexedMoments<Dim,3>::size();
  }
  else if(a_order == 4)
  {
    retval = IndexedMoments<Dim,4>::size();
  }
  else
  {
    amrex::Error("need to add more cases to getIndMomSize");
  }
  return retval;
}

///
/**
   Moments are centered at the center of the cell.
   For each of these moments I shift them to the lowest 
   corner of the cell, where I know what the bounds of the integrals 
   is (lower bound always zero, upper bound = dx^d dx^px dx^py dx^pz
   If the shifted moment is out of bounds, I bound it.
   The tolerance is about verbosity.    If the moment is outside the tolerance
   then it gets included into a_bogusPowers.
   EBMoments do not get checked for maxvals.
**/
template<int Dim, int ORDER>
void
checkMoments(IndexedMoments<Dim, ORDER> & a_moments,
             Vector<IndexTM<int, Dim> > & a_bogusPowers,
             const Real                 & a_dx,
             const Real                 & a_tolerance,
             const bool                 & a_ebMoment,
             const bool                 & a_bindMoments)
{
  MomentIterator<Dim, ORDER> momit;
  a_bogusPowers.resize(0);
  Real cellvol = 1;
  for(int idir = 0; idir < Dim; idir++)
  {
    cellvol *= a_dx;
  }

  //shift moments to lowest corner 
  IndexTM<Real, Dim> shiftvec;
  shiftvec.setAll(0.5*a_dx);
  a_moments.shift(shiftvec);

  //because they are all shifted, to the lowest 
  //min is always zero.
  const Real zerval = 0;
  for(momit.reset(); momit.ok(); ++momit)
  {
    const IndexTM<int, Dim>& p = momit();
    //this makes max = dx^p
    Real maxval = POW<Dim>(a_dx, p);
    //because this is a vol integral...
    maxval *= cellvol;
    //this is a reference so it will change the input
    Real& momval = a_moments[p];
    if(momval < zerval)
    {
      if((zerval - momval) > a_tolerance*maxval)
      {
        a_bogusPowers.push_back(p);
        //              MayDay::Error("minval break low");
      }
      if(a_bindMoments)
      {
        momval = zerval;
      }
    }
    else if((momval > maxval) && (!a_ebMoment))
    {
      if((momval - maxval) > a_tolerance*maxval)
      {
        a_bogusPowers.push_back(p);
        //              MayDay::Error("maxval break");
      }
      if(a_bindMoments)
      {
        momval = maxval;
      }
    }
  }

  //shift moments back to highest corner 
  shiftvec.setAll(-a_dx);
  a_moments.shift(shiftvec);

  //because the moments are shifted to the highest corner
  //even powers will always be positive, odd ones will always be negtive
  for(momit.reset(); momit.ok(); ++momit)
  {
    const IndexTM<int, Dim>& p = momit();
    int psum = p.sum();
    bool pEven = (psum % 2 == 0);
    Real& momval = a_moments[p];
    //this makes max = dx^p
    Real maxval = POW<Dim>(a_dx, p);
    //because this is a vol integral...
    maxval *= cellvol;
    if(pEven && (momval < zerval) && (!a_ebMoment))
    {
      //p has an  sum and therefore the moment should be positive
      if(momval < -a_tolerance*maxval )
      {
        a_bogusPowers.push_back(p);
        //              MayDay::Error("high corner break even");
      }
      if(a_bindMoments)
      {
        momval = 0;
      }
    }
    else if(!pEven && (momval > zerval) && (!a_ebMoment))
    {
      //p has an odd sum and therefore the moment should be negative
      if(momval > a_tolerance*maxval)
      {
        a_bogusPowers.push_back(p);
        //              MayDay::Error("high corner break odd");
      }
      if(a_bindMoments)
      {
        momval = 0;
      }
    }
  }

  //shift moments back to cell center
  shiftvec.setAll(0.5*a_dx);
  a_moments.shift(shiftvec);
}


///
/**
   return true if all of a_index >= 0, false otherwise
*/
template <int Dim>
bool allPositive(const IndexTM<int, Dim>& a_index)
{
  bool allpos = true;
  for(int idir = 0; idir < Dim; idir++)
  {
    if(a_index[idir] < 0)
    {
      allpos = false;
    }
  }
  return allpos;
}


#endif
