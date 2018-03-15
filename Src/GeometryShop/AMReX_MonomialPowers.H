#include "AMReX_IndexTM.H"
#include "AMReX_Vector.H"
#ifndef _MonomialPowers_H_
#define _MonomialPowers_H_

#define CH_IM_MAX_POWER  8
using amrex::Vector;
static const int s_mpowers_max_sizes[3][CH_IM_MAX_POWER+1] =
{
  {1, 2,  3,  4,  5,  6,  7,   8,   9}, // D=1
  {1, 3,  6, 10, 15, 21, 28,  36,  45}, // D=2
  {1, 4, 10, 20, 35, 56, 84, 120, 165}  // D=3
};

template <int Dim, int P>
class MonomialPowers
{
public:

  ///number of reals in the vector
  static int size()
    {
      if(!s_staticsSet) 
      {
        setStatics();
      }
      return s_size;
    }

  ///monomial powers 
  static const Vector<IndexTM<int,Dim> >& getMonomialPowers()
    {
      if(!s_staticsSet) 
      {
        setStatics();
      }
      return s_multiIndicies;
    }


protected:

  ///
  static void setStatics()
    {
      BL_ASSERT(Dim <= BL_SPACEDIM);
      BL_ASSERT(Dim > 0);
      BL_ASSERT(P <  CH_IM_MAX_POWER);
      BL_ASSERT(P >= 0);

      s_size = s_mpowers_max_sizes[Dim-1][P];

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

      s_staticsSet = true;
    }
  
  ///
  static bool s_staticsSet;

  ///
  static  int s_size;

  ///
  static Vector<IndexTM<int,Dim> > s_multiIndicies;

  
private:

    
}; 

template <int Dim, int P> bool                       MonomialPowers<Dim, P>::s_staticsSet     = false;
template <int Dim, int P> int                        MonomialPowers<Dim, P>::s_size           = -1;
template <int Dim, int P> Vector<IndexTM<int,Dim> >  MonomialPowers<Dim, P>::s_multiIndicies  = Vector<IndexTM<int,Dim> >();

#endif
