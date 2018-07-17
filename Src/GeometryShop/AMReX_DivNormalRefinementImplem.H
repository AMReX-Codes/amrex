
#ifndef AMREX_DIVNORMALREFINEMENTIMPLEM_H_
#define AMREX_DIVNORMALREFINEMENTIMPLEM_H_

template <int dim> DivNormalRefinement<dim>::DivNormalRefinement()
{
}

template <int dim> DivNormalRefinement<dim>::DivNormalRefinement(const Real & a_maxChangeThreshold,
                                                                 const int& a_maxNumberOfRefinements)
{
  setMaxChangeThreshold(a_maxChangeThreshold);
  setMaxNumberOfRefinements(a_maxNumberOfRefinements);
}

template <int dim> DivNormalRefinement<dim>::~DivNormalRefinement()
{
}

template <int dim> bool DivNormalRefinement<dim>::doRefine(IndexTM<int,dim>          & a_refineInDir,
                                                           const CutCellMoments<dim> & a_ccm,
                                                           const int                 & a_numberOfRefinements)
{
  bool retval = false;
  a_refineInDir = IndexTM<int,dim>::Zero;

  if (a_numberOfRefinements < m_maxNumberOfRefinements)
  {
    Real changeInNormal = approximateDivNormal(a_ccm);

    //check whether normal equals the zero vector or whether constraints were active
    if (changeInNormal > m_maxChangeThreshold)
    {
      retval        = true;
      a_refineInDir = IndexTM<int,dim>::Unit;
    }
  }

  return retval;
}

template <int dim> Real DivNormalRefinement<dim>::approximateDivNormal(const CutCellMoments<dim> & a_ccm)
{

  typedef IndexTM<Real,dim> RvDim;
  typedef IndexTM<int,dim>                      EdgeIndex;
  typedef map<EdgeIndex,Real,LexLT<EdgeIndex> > EdgeIntersections;

  Real changeInNormal = 0.0;
  const IFData<dim>& edgeData = a_ccm.m_IFData;

  //normal at the center of the Taylor Expansion, which always exists.
  const RvDim& normal = edgeData.m_normalDerivatives.find(IndexTM<int,dim>::Zero)->second;

  //iterate through intersection points
  for (typename EdgeIntersections::const_iterator it = edgeData.m_intersections.begin();it != edgeData.m_intersections.end(); ++it)
  {
    const EdgeIndex& edgeIndex = it->first;
    const Real&      intercept = it->second;

    int varyDir = edgeIndex[0];

    RvDim intersectPt = RvDim::Zero;
    intersectPt[varyDir] = intercept;

    for (int i = 1; i < dim; i++)
    {
      int curDir;
      int loHi;

      if (i <= varyDir)
      {
        curDir = i-1;
      }
      else
      {
        curDir = i;
      }

      loHi = edgeIndex[i];
      intersectPt[curDir]  = loHi - 0.5;
      intersectPt[curDir] *= edgeData.m_globalCoord.m_dx[curDir];
    }

    //normal at intersection point is the zeroth derivative of the normal
    NormalDerivative<dim> normalDerivative;
    RvDim intersectPtNormal;
    for (int idir = 0; idir < dim; idir++)
    {
      intersectPtNormal[idir] = normalDerivative.evaluate(IndexTM<int,dim>::Zero,
                                                          idir,
                                                          edgeData.m_globalCoord.convert(intersectPt,edgeData.m_cellCenterCoord),
                                                          edgeData.m_function);
    }

    // find the angle between the normal at intersection point the normal at the average of interesection points
    Real dotProduct = 0.0;
    for (int idir = 0; idir < dim; idir++)
    {
      dotProduct += intersectPtNormal[idir]*normal[idir];
    }

    if (std::abs(dotProduct - 1.0) > changeInNormal)
    {
      changeInNormal = std::abs(dotProduct -1.0);
    }
  }
  return changeInNormal;
}

template <int dim> void DivNormalRefinement<dim>::setMaxChangeThreshold(const Real & a_maxChangeThreshold)
{
  if (a_maxChangeThreshold < 0)
  {
    amrex::Abort("DivNormalRefinement<dim>::setMaxChangeThreshold - maxChangeThreshold must be >= 0");
  }

  m_maxChangeThreshold = a_maxChangeThreshold;
}

template <int dim> Real DivNormalRefinement<dim>::getMaxChangeThreshold()
{
  return m_maxChangeThreshold;
}

template <int dim> void DivNormalRefinement<dim>::setMaxNumberOfRefinements(const int & a_maxNumberOfRefinements)
{
  if (a_maxNumberOfRefinements < 0)
  {
    amrex::Abort("DivNormalRefinement<dim>::setMaxNumberOfRefinements - maxNumberOfRefinents must be >= 0");
  }

  m_maxNumberOfRefinements = a_maxNumberOfRefinements;
}

template <int dim> int DivNormalRefinement<dim>::getMaxNumberOfRefinements()
{
  return m_maxNumberOfRefinements;
}

#endif
