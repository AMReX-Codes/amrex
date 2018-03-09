
#ifndef _INDEXEDMOMENTS_H_
#define _INDEXEDMOMENTS_H_

#include "AMReX_IndexTM.H"
#include "AMReX_Vector.H"
#include "AMReX_SPACE.H"
#include "AMReX_CH_EBIS_ORDER.H"


namespace amrex
{
  using std::map;

#define CH_IM_MAX_POWER  8

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
    IndexedMoments();


    //! Destructor
    ~IndexedMoments() {};


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
    IndexedMoments<Dim, P>&  operator+=(const IndexedMoments<Dim, P>& increment);

    ///multiply each component by constant
    IndexedMoments<Dim, P>&  operator*=(const Real& a_factor);


    ///number of reals in the vector
    static int size()
      {
        if(!s_staticsSet) 
        {
          setStatics();
        }
        return s_size;
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
    void setRegular(const Real a_dx);

    ///monomial powers 
    static const Vector<IndexTM<int,Dim> >& getMonomialPowers()
      {
        if(!s_staticsSet) 
        {
          setStatics();
        }
        return s_multiIndicies;
      }

    ///
    bool isRegular() const
      {
        return (m_isRegular==1);
      }

    /// for use with irregnode
    IndexedMoments<Dim, P>& operator=(const map<IndexTM<int, Dim>,  Real>& a_mapin);

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
    void shift(const IndexTM<Real, Dim>& a_distance);
  
    ///
    void setToZero()
      {
        for(int ivec = 0; ivec < s_size; ivec++)
        {
          m_moms[ivec] = 0.;
        }
      }

    /// Calculate what linear index this multi-index is
    static int indexOf(IndexTM<int,Dim> a_index);

    ///
    static IndexTM<int,Dim> getIndex(const int& a_linearIndex)
      {
        return s_multiIndicies[a_linearIndex];
      }

    /// outputs the current state to pout() (a la parstream.H)
    void spout() const;

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
                                const IndexedMoments<Dim, P> & a_CB);



    ///divides each entry by p!
    void divideByFactorial();

    ///multiply each entry by p!
    void multiplyByFactorial();

  protected:

    ///
    Real 
    getMoment(const IndexTM<int, Dim>        & a_mono,
              const map<IndexTM<int, Dim>,  Real>& a_mapin) const;

    ///
    static void setStatics();
  
    ///
    static bool s_staticsSet;

    ///
    static  int s_size;

    ///
    static Vector<IndexTM<int,Dim> > s_multiIndicies;

    ///
    static void setMultiIndicies();

    ///
    static void setSize();
  
  private:

    // Indicator that we contain only "full" moments
    int m_isRegular; 

    // Our moments to store    
    Vector<Real> m_moms;

    static const int s_max_sizes[][CH_IM_MAX_POWER+1];
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
      MayDay::Error("need to add more cases to getIndMomLinearIndex");
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
      MayDay::Error("need to add more cases to getIndMomMultiIndex");
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
      MayDay::Error("need to add more cases to getIndMomSize");
    }
    return retval;
  }

///to see if all powers of p are even
  template<int Dim>
  bool allEven(const IndexTM<int, Dim>& a_p)
  {
    bool retval = true;
    for(int idir = 0; idir < Dim; idir++)
    {
      if(a_p[idir]%2 != 0)
      {
        retval = false;
        break;
      }
    }
    return retval;
  }

/// computes x^p
  template<int Dim>
  Real 
  POW(const Real& a_x,  const IndexTM<int, Dim> a_p)
  {
    Real retval = 1;
    for(int idir = 0; idir < Dim; idir++)
    {
      for(int iexp = 0; iexp < a_p[idir]; iexp++)
      {
        retval *= a_x;
      }
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
               const bool                 & a_bindMoments);

///
/**
   return true if all of a_index >= 0, false otherwise
*/
  template <int Dim>
  bool allPositive(const IndexTM<int, Dim>& a_index); 

/**/
///template specializations for debugger
/**
   template< >
   void
   checkMoments(IndexedMoments<SpaceDim, CH_EBIS_ORDER> & a_moments,
   Vector<IndexTM<int, SpaceDim> > & a_bogusPowers,
   const Real                 & a_dx,
   const Real                 & a_tolerance,
   const bool                 & a_ebMoment,
   const bool                 & a_bindMoments);


   template< >
   void
   checkMoments(IndexedMoments<SpaceDim-1, CH_EBIS_ORDER> & a_moments,
   Vector<IndexTM<int, SpaceDim-1> > & a_bogusPowers,
   const Real                 & a_dx,
   const Real                 & a_tolerance,
   const bool                 & a_ebMoment,
   const bool                 & a_bindMoments);
**/

}

#include "AMReX_IndexedMomentsImplem.H"

#endif
