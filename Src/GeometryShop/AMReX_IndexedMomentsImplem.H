#ifndef _IndexedMomentsImplem_H_
#define _IndexedMomentsImplem_H_
#include "MomentIterator.H"
#include "EB_TYPEDEFS.H"
#include "IndexTM.H"
#include "Factorial.H"
#include "CH_Timer.H"
#include "parstream.H"
#include <map>
#include <utility>
#include "NamespaceHeader.H"

template <int Dim, int P> bool                      IndexedMoments<Dim, P>::s_staticsSet     = false;
template <int Dim, int P> Vector<IndexTM<int,Dim> > IndexedMoments<Dim, P>::s_multiIndicies  = Vector<IndexTM<int, Dim> >();
template <int Dim, int P> int                       IndexedMoments<Dim, P>::s_size           = 0;
using std::endl;

///
/**
   Moments are centered at the center of the cell.
   For each of these moments I shift them to the lowest 
   corner of the cell, where I know what the bounds of the integrals 
   is (lower bound always zero, upper bound = dx^d dx^px dx^py dx^pz
   If the shifted moment is out of bounds, I bound it throw a message.
   tolerance is about verbosity.  
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

template <int Dim, int P>
const int IndexedMoments<Dim, P>::s_max_sizes[/* D */][CH_IM_MAX_POWER+1] = 
  {
    {1, 2,  3,  4,  5,  6,  7,   8,   9}, // D=1
    {1, 3,  6, 10, 15, 21, 28,  36,  45}, // D=2
    {1, 4, 10, 20, 35, 56, 84, 120, 165}  // D=3
  }; 
using std::endl; 
/////
template <const int Dim, const int P>
Real
IndexedMoments<Dim, P>:: 
getMoment(const     IndexTM<int, Dim>        & a_mono,
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
      std::cout << endl;
      MayDay::Abort("getMoments in IndexedMoments: moment not found in map");
    }
  return moment;
}

///multiply each component by factor
template <int Dim, int P>
IndexedMoments<Dim, P>&
IndexedMoments<Dim, P>:: 
operator*=(const Real& a_factor)
{
  for(MomentIterator<Dim,P> it; it.ok(); ++it)
    {
      (*this)[it()] *= a_factor;
    }
  return *this;
}
///add compoenentwise
template <int Dim, int P>
IndexedMoments<Dim, P>&
IndexedMoments<Dim, P>:: 
operator+=(const IndexedMoments<Dim, P>& a_input)
{
  for(MomentIterator<Dim,P> it; it.ok(); ++it)
    {
      (*this)[it()] += a_input[it()];
    }
  return *this;
}
///shift moment by the input distance.
template <int Dim, int P>
void
IndexedMoments<Dim, P>:: 
shift(const IndexTM<Real, Dim>& a_distance)
{
  CH_TIME("IndexedMoment::shift");
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
                  xpow *=  POW(x[idir], p[idir]);
                }
              
              moment += pCk<Dim>(r, q) * xpow * m0;
            
            } //end if (q <= r)
        }//end loop over inner moments
      retval[r] = moment;
    }
  *this = retval;
}

/*******/
/// for use with irregnode
template <int Dim, int P>
IndexedMoments<Dim, P>& 
IndexedMoments<Dim, P>:: 
operator=(const map<IndexTM<int, Dim>,  Real>& a_mapin)
{
  
  for(MomentIterator<Dim, P> iter; iter.ok(); ++iter)
    {
      IndexTM<int, Dim> mono = iter();
      Real moment= getMoment(iter(), a_mapin);
      (*this)[iter()] = moment;
    }
  return *this;
}

template <int Dim, int P>
void
IndexedMoments<Dim, P>:: 
spout() const
{
  pout() << "index \t  moment  " << endl;
  for(MomentIterator<Dim, P> iter; iter.ok(); ++iter)
    {
      pout() <<  iter() << "\t" << (*this)[iter()] << endl;
    }
}

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


template <int Dim, int P>
int
IndexedMoments<Dim, P>::
indexOf(IndexTM<int,Dim> a_index)
{
  CH_assert(a_index.sum() <= P);
  CH_assert(allPositive(a_index));

  // Formula for the layout index with 0th dim first, then 1st, etc.
  int index = a_index[0];
  int rem = P;
  for (int d=Dim-1; d > 0; --d)
    {
      index += s_max_sizes[d][rem] - s_max_sizes[d][rem-a_index[d]];
      rem -= a_index[d];
    }

  return index;
}
/*******/
template <int Dim, int P>
void
IndexedMoments<Dim, P>::
setStatics() 
{
  setSize();
  setMultiIndicies();
  s_staticsSet= true;
}

/*******/
template <int Dim, int P>
IndexedMoments<Dim, P>::
IndexedMoments() 
{
  if(!s_staticsSet)
    {
      setStatics();
    }

  m_moms.resize(s_size);
  m_isRegular = 0;
}

/*******/
template <int Dim, int P>
void
IndexedMoments<Dim, P>::
setRegular(const Real a_dx)
{
  CH_assert(a_dx > 0);
  m_isRegular = 1;

  // Get the set of multi-indexes in the correct order
  const Vector<IndexTM<int,Dim> >& mix = IndexedMoments<Dim, P>::getMonomialPowers();

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
          moment = POW(0.5 * a_dx, index.sum()) * POW(a_dx, Dim);
          for (int d=0; d < Dim; ++d)
            moment /= (Real) (index[d]+1);
        }
      m_moms[ix] = moment;

      // pout() << ix << ": " << index << " = " << moment << endl;
    }

}
/*******/
template <int Dim, int P>
void
IndexedMoments<Dim, P>::
setSize()
{
  CH_assert(Dim <= SpaceDim);
  CH_assert(Dim > 0);
  CH_assert(P <  CH_IM_MAX_POWER);
  CH_assert(P >= 0);
  

  s_size = s_max_sizes[Dim-1][P];
}
/*******/
template <int Dim, int P>
void
IndexedMoments<Dim, P>::
setMultiIndicies()
{
  s_multiIndicies.resize(s_size);

  IndexTM<int,Dim> index = IndexTM<int,Dim>::Zero;
  for (int ix=0; ix < s_size; ++ix)
    {
      // If the sum is too large, shift extras to the right
      for (int d=0; (index.sum() > P) && (d < Dim-1); ++d)
        {
          index[d] = 0;
          ++index[d+1];
        }
      s_multiIndicies[ix] = index;
      ++index[0];
    }
}
/****/
template <int Dim, int P>
void 
IndexedMoments<Dim, P>::
setToTruncatedMultiply(const IndexedMoments<Dim, P> & a_CA,
                       const IndexedMoments<Dim, P> & a_CB)
{
  this->setToZero();
  //here we are computing sum(m_p) = integral(ni(x-xbar)^p1/p! (x-xbar)^p2)
  // = sum_p2 ((n^p2)/p2!) m_p+p12
  for(MomentIterator<Dim, P> momitOuter; momitOuter.ok(); ++momitOuter)
    {
      //this is where the answer goes 
      IndexTM<int,SpaceDim>  po = momitOuter();

      for(MomentIterator<Dim, P> momitInner; momitInner.ok(); ++momitInner)
        {
          IndexTM<int,SpaceDim> pi = momitInner();
          IndexTM<int,SpaceDim> psum = po + pi;
          if(psum.sum() <= P)
            {
              Real incr = a_CA[pi]*a_CB[po];
              (*this)[psum] += incr;
            }
        }
    }

}
/****/
template <int Dim, int P>
void 
IndexedMoments<Dim, P>::
divideByFactorial()
{
  for(MomentIterator<Dim, P> momit; momit.ok(); ++momit)
    {
      IndexTM<int,SpaceDim>  p = momit();
      Real factorial = pfactorial(p);
      (*this)[p] /= factorial;
    }
}
/****/
template <int Dim, int P>
void 
IndexedMoments<Dim, P>::
multiplyByFactorial()
{
  for(MomentIterator<Dim, P> momit; momit.ok(); ++momit)
    {
      IndexTM<int,SpaceDim>  p = momit();
      Real factorial = pfactorial(p);
      (*this)[p] *= factorial;
    }
}
/*******/
#include "NamespaceFooter.H"
#endif
