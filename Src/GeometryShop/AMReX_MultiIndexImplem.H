
#ifndef _MULTIINDEXIMPLEM_H_
#define _MULTIINDEXIMPLEM_H_


// Generate and return a vector of all multi-indices whose magnitude (the sum
// of their entries) equals "a_magnitude"
template <int dim> void generateMultiIndices(Vector<IndexTM<int,dim> > & a_indices,
                                             const int                 & a_magnitude)
{
  if (a_magnitude >= 0)
  {
    IndexTM<int,dim> monomial;

    for (int idir = 0; idir < dim; ++ idir)
    {
      monomial[idir] = 0;
    }

    while (true)
    {
      for (int j = 1; j < dim-1; ++j)
      {
        int t = a_magnitude;
        for (int k = j+1; k < dim; ++k)
        {
          t -= monomial[k];
        }

        if (monomial[j] > t)
        {
          monomial[j+1] += 1;
          monomial[j] = 0;
        }
        else
        {
          break;
        }
      }

      if (monomial[dim-1] > a_magnitude)
      {
        break;
      }

      monomial[0] = a_magnitude;

      for (int j = 1; j < dim; ++j)
      {
        monomial[0] -= monomial[j];
      }

      a_indices.push_back(monomial);

      monomial[1] += 1;
    }
  }
}


#endif
