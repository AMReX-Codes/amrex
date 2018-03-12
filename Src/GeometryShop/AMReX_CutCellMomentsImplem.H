
#ifndef _CUTCELLMOMENTSIMPLEM_H_
#define _CUTCELLMOMENTSIMPLEM_H_

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <string>


#include "AMReX_LSProblem.H"

namespace amrex
{
  using std::string;
// Null constructor
  template <int dim> CutCellMoments<dim>::CutCellMoments()
  {
  }

// Copy constructor
  template <int dim> CutCellMoments<dim>::CutCellMoments(const CutCellMoments<dim>& a_cutCellMoments)
    :m_moments         (a_cutCellMoments.m_moments),
     m_EBmoments       (a_cutCellMoments.m_EBmoments),
     m_bdCutCellMoments(a_cutCellMoments.m_bdCutCellMoments),
     m_IFData          (a_cutCellMoments.m_IFData),
     m_bdCCOn          (a_cutCellMoments.m_bdCCOn),
     m_residual        (a_cutCellMoments.m_residual),
     m_numActiveBounds (a_cutCellMoments.m_numActiveBounds),
     m_badNormal       (a_cutCellMoments.m_badNormal)
  {
  }

// Use this for initializing
  template <int dim> CutCellMoments<dim>::CutCellMoments(const IFData<dim>& a_info)
    :m_IFData(a_info),
     m_numActiveBounds(0)
  {
    m_bdCCOn    = false;
    m_badNormal = false;

    for (int hilo = 0; hilo < 2; ++hilo)
    {
      for (int idir = 0; idir < dim; ++idir)
      {
        //identifier for which boundary cutCellMoment
        Iv2 bdId;
        bdId[BDID_DIR] = idir;
        bdId[BDID_HILO] = hilo;

        IFData<dim-1> reducedInfo(a_info,a_info.m_maxOrder+1,idir,hilo);
        CutCellMoments<dim-1>bdCutCellMoments(reducedInfo);

        m_bdCutCellMoments[bdId] = bdCutCellMoments;


        // Notice whether reducedInfo coincides with the interface
        if (reducedInfo.m_allVerticesOn)
        {
          m_bdCCOn = true;
        }

        // Notice whether reducedInfo.m_normal = zero vector
        if (reducedInfo.m_badNormal)
        {
          m_badNormal = true;
        }

      }
    }
  }

// Destructor
  template <int dim> CutCellMoments<dim>::~CutCellMoments()
  {
  }

  template <int dim> const CutCellMoments<dim - 1> CutCellMoments<dim>::getBdCutCellMoments(const Iv2 & a_bdId) const
  {
    typename BdCutCellMoments::const_iterator it = m_bdCutCellMoments.find(a_bdId);

    if (it == m_bdCutCellMoments.end())
    {
      amrex::Abort("Can't find this bdId in m_bdCutCellMoments");
    }

    return it->second;
  }

// This function is called to fill the boundary moments from the higher
// dimension: it fills the moments in a_coarseCutCell using either the
// moments in the current cutcell map or the full cell/covered cell answers
// depending on the value of a_IFData.m_allVerticesIn/Out
  template <int dim> void CutCellMoments<dim>::addBdMoments(CutCellMoments<dim>     & a_coarseBdCutCell,
                                                            const IFData<dim+1>     & a_IFData,
                                                            const int               & a_degreePmax,
                                                            const bool              & a_useConstraints,
                                                            const IndexTM<Real,dim> & a_refinedCenterDelta,
                                                            const IndexTM<int,dim>  & a_hilo)
  {
    typedef map<IvDim,int > PthMomentLoc;

    if (a_IFData.m_allVerticesIn)
    {
      // The boundary map doesn't exist, we need to compute the full cell
      // quadrature for all monomial in order to use the changecoordinates
      // function.  Generate a full cell quadrature map.
      PthMoment fullCellMap;

      for (int iDegree = a_degreePmax; iDegree > 0; iDegree--)
      {
        // Creating a lsproblem generates the monomial map that is
        // needed to fill the moment maps
        LSProblem<dim> lsp(iDegree,a_useConstraints);
        const PthMomentLoc& momMap = lsp.getMonomialLocMapDegreePLess1();

        for (typename PthMomentLoc::const_iterator it = momMap.begin(); it != momMap.end(); ++it)
        {
          //        pout()<<"mono = "<<it->first<<endl;
          fullCellMap[it->first] = fullCellQuadrature(it->first,m_IFData.m_cellCenterCoord);
        }
      }

      for (int iDegree = a_degreePmax; iDegree > 0; iDegree--)
      {
        LSProblem<dim> lsp(iDegree,a_useConstraints);

        // Change the coordinates of the refined moments using the full
        // cell quadrature map.  Add up the obtained moment to the
        // coarse map
        const PthMomentLoc& momMap = lsp.getMonomialLocMapDegreePLess1();

        for (typename PthMomentLoc::const_iterator it = momMap.begin(); it != momMap.end(); ++it)
        {
          Real bdMoment = getBdMoment(it->first,a_IFData,a_refinedCenterDelta,fullCellMap);
          a_coarseBdCutCell.m_moments[it->first] += bdMoment;
        }
      }
    }
    else
    {
      for (int iDegree = a_degreePmax; iDegree > 0; iDegree--)
      {
        LSProblem<dim> lsp(iDegree,a_useConstraints);

        // Get the moments from the refined cell map + change coordinates
        // function if the refined cell is irregular or they are 0 as all
        // vertices are out
        const PthMomentLoc& momMap = lsp.getMonomialLocMapDegreePLess1();

        for (typename PthMomentLoc::const_iterator it = momMap.begin(); it != momMap.end(); ++it)
        {
          Real bdMoment = getBdMoment(it->first,a_IFData,a_refinedCenterDelta);
          a_coarseBdCutCell.m_moments[it->first] += bdMoment;
        }
      }
    }

    for (int iDegree = a_degreePmax; iDegree >= 0; iDegree--)
    {
      LSProblem<dim> lsp(iDegree,a_useConstraints);

      // EB moments are either found in the EB moment map or equals to 0
      // (allIn or allOut)
      const PthMomentLoc& momMapP = lsp.getMonomialLocMapDegreeP();

      for (typename PthMomentLoc::const_iterator it = momMapP.begin(); it != momMapP.end(); ++it)
      {
        Real bdEBMoment = getBdEBMoment(it->first,a_IFData,a_refinedCenterDelta);
        a_coarseBdCutCell.m_EBmoments[it->first] += bdEBMoment;
      }
    }

    IndexTM<int,2> bdId;

    for (int idir = 0; idir < dim; idir++)
    {
      bdId[0] = idir;

      // The refined cell center in dim-1 is needed to change the coordinate
      // system of the refined moments to the coordinate system where they are
      // needed.
      IndexTM<Real,dim-1> localRefinedCellCenter;
      IndexTM<int,dim-1> localHilo;

      for (int jdir = 0; jdir < dim; jdir++)
      {
        if (jdir < idir)
        {
          localRefinedCellCenter[jdir] = a_refinedCenterDelta[jdir];
          localHilo[jdir] = a_hilo[jdir];
        }
        else if (jdir > idir)
        {
          localRefinedCellCenter[jdir-1] = a_refinedCenterDelta[jdir];
          localHilo[jdir-1] = a_hilo[jdir];
        }
      }

      if (a_hilo[idir] != LARGEINTVAL)
      {
        bdId[1] = a_hilo[idir];
        // Drop one dimension to add up the boundary moments. Dropping one
        // dimension allows an easier implementation as 1D is a particular
        // case
        (m_bdCutCellMoments[bdId]).addBdMoments(a_coarseBdCutCell.m_bdCutCellMoments[bdId],m_IFData,a_degreePmax+1,a_useConstraints,localRefinedCellCenter,localHilo);
      }
      else
      {
        for (int side = 0; side < 2; side++)
        {
          bdId[1]=side;
          (m_bdCutCellMoments[bdId]).addBdMoments(a_coarseBdCutCell.m_bdCutCellMoments[bdId],m_IFData,a_degreePmax,a_useConstraints,localRefinedCellCenter,localHilo);
        }
      }
    }
  }

  template <int dim> Real CutCellMoments<dim>::fullCellQuadrature(const IndexTM<int,dim>      & a_mono,
                                                                  const CoordinateSystem<dim> & a_coord)
  {
    Real moment = 1.0;

    for (int idir = 0; idir < dim; ++idir)
    {
      Real loPt = a_coord.convertDir(-0.5*m_IFData.m_cellCenterCoord.m_dx[idir],
                                     m_IFData.m_cellCenterCoord,
                                     idir);

      Real hiPt = a_coord.convertDir(0.5*m_IFData.m_cellCenterCoord.m_dx[idir],
                                     m_IFData.m_cellCenterCoord,
                                     idir);

      Real idirInt = pow(hiPt,a_mono[idir] + 1) - pow(loPt,a_mono[idir] + 1);

      idirInt /= (a_mono[idir] + 1);
      moment *= idirInt;
    }

    return moment;
  }

  template <int dim> Real CutCellMoments<dim>::changeMomentCoordinates(PthMoment        & a_refinedMomentMap,
                                                                       const IndexTM<int,dim>  & a_monomial,
                                                                       const IndexTM<Real,dim> & a_refinedCenterDelta)
  {
    // This function outputs the moment in the coordinate system where the
    // origin is the current cell center, given the moment map of the refined
    // moments computed in the coordinate system where the origin is the
    // refined cell center it uses a Taylor development around 0 of the function
    // (x+xc)^p1*(y+yc)^p2*...
    Real moment = 0.0;

    // Degree is the degree of the monomial
    int degree = 0;

    for (int i = 0; i < dim; i++)
    {
      degree += a_monomial[i];
    }

    // Loop over the possible values of r in the development of the Taylor series
    // The Taylor series has a finite expansion as the function's derivatives
    // are zero for r>degree
    for (int r = 0; r <= degree; r++)
    {
      if (r >= 0)
      {
        // Generate all the possible monomials of degree r and add up the moment
        IndexTM<int,dim> derivative;

        for (int idir = 0; idir < dim; ++idir)
        {
          derivative[idir] = 0;
        }

        while (true)
        {
          for (int j = 1; j < dim-1; ++j)
          {
            int t = r;
            for (int k = j+1; k < dim; ++k)
            {
              t -= derivative[k];
            }

            if (derivative[j] >  t)
            {
              derivative[j+1] += 1;
              derivative[j  ]  = 0;
            }
            else
            {
              break;
            }
          }

          if (derivative[dim-1] > r)
          {
            break;
          }

          derivative[0] = r;

          for (int j = 1; j < dim; ++j)
          {
            derivative[0] -= derivative[j];
          }

          // Add up the relevant term of the refined moment map to the map
          Real coeff = 1.0;
          Real factorial = 1.0;

          for (int idir = 0; idir < dim; idir++)
          {
            for (int j = 0; j < derivative[idir]; j++)
            {
              coeff *= a_monomial[idir] - j;
              factorial *= j+1;
            }

            // This is to prevent the computation of 0*0^-1 which result
            // in a nan
            if (a_monomial[idir]-derivative[idir] < 0)
            {
              if (coeff != 0)
              {
                amrex::Abort("Coeff should be zero when the derivative has higher coefficient than the monomial");
              }
            }
            else
            {
              coeff *= pow(a_refinedCenterDelta[idir],a_monomial[idir]-derivative[idir]);
            }
          }

          //pout()<<a_refinedMomentMap[derivative]<<endl;
          moment += coeff * a_refinedMomentMap[derivative] / factorial;

          // Increments to get to the next derivative
          derivative[1] += 1;
        }
      }
    }

    return moment;
  }

  template <int dim> void CutCellMoments<dim>::changeMomentCoordinatesToCellCenter()
  {
    // Move moments from parent coord to cell center coord
    IndexTM<Real,dim> delta = m_IFData.m_cellCenterCoord.m_origin;
    delta                  -= m_IFData.m_parentCoord    .m_origin;

    PthMoment copyMoments = m_moments;
    for (typename PthMoment::const_iterator it = copyMoments.begin();
         it != copyMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_moments[mono] = changeMomentCoordinates(copyMoments, mono, delta);
    }

    PthMoment copyEBMoments = m_EBmoments;
    for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_EBmoments[mono] = changeMomentCoordinates(copyEBMoments, mono, delta);
    }
  }

  template <int dim> void CutCellMoments<dim>::changeMomentCoordinatesToParentCenter()
  {
    // Move moments from cell center coord to parent coord
    IndexTM<Real,dim> delta = m_IFData.m_parentCoord    .m_origin;
    delta                  -= m_IFData.m_cellCenterCoord.m_origin;

    PthMoment copyMoments = m_moments;
    for (typename PthMoment::const_iterator it = copyMoments.begin();
         it != copyMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_moments[mono] = changeMomentCoordinates(copyMoments, mono, delta);
    }

    PthMoment copyEBMoments = m_EBmoments;
    for (typename PthMoment::const_iterator it = copyEBMoments.begin(); it != copyEBMoments.end(); ++it)
    {
      IvDim mono = it->first;
      m_EBmoments[mono] = changeMomentCoordinates(copyEBMoments, mono, delta);
    }
  }

  template <int dim> void CutCellMoments<dim>::initialize(CutCellMoments<dim> & a_refinedCutCell)
  {
    // This initialize the maps of the current cut cell moment to zero by
    // iterating through the refined maps
    initializeMap(m_EBmoments,a_refinedCutCell.m_EBmoments);
    initializeMap(m_moments,a_refinedCutCell.m_moments);

    // Initialize maps on the boundary
    IndexTM<int,2> bdId;

    for (int idir = 0; idir < dim; idir++)
    {
      bdId[0] = idir;

      for (int hilo = 0; hilo < 2; hilo++)
      {
        bdId[1] = hilo;
        CutCellMoments<dim-1> refinedBdCutCell = a_refinedCutCell.getBdCutCellMoments(bdId);

        m_bdCutCellMoments[bdId].initialize(refinedBdCutCell);
      }
    }
  }

  template <int dim> void CutCellMoments<dim>::initializeMap(PthMomentLesserDimension & a_map1,
                                                             PthMomentLesserDimension & a_map2)
  {
    for (typename PthMomentLesserDimension::const_iterator it = a_map2.begin(); it != a_map2.end(); ++it)
    {
      a_map1[it->first] = 0.0;
    }
  }

  template <int dim> void CutCellMoments<dim>::initializeMap(PthMoment & a_map1,
                                                             PthMoment & a_map2)
  {
    for (typename PthMoment::const_iterator it = a_map2.begin(); it != a_map2.end(); ++it)
    {
      a_map1[it->first] = 0.0;
    }
  }


  template <int dim> Real CutCellMoments<dim>::getBdMoment(const IvDim             & a_mono,
                                                           const IFData<dim+1>     & a_IFData,
                                                           const IndexTM<Real,dim> & a_refinedCenterDelta,
                                                           PthMoment                 a_fullCellMap)
  {
    // This computes the moments on the fly in the case the refined cell (in
    // dim+1) was all in or all out (in that case the boundary moments were
    // never computed) or looks for them in the moment map
    Real moment = LARGEREALVAL;

    if (a_IFData.m_allVerticesOut)
    {
      moment = 0.0;
    }
    else if (a_IFData.m_allVerticesIn)
    {
      moment = changeMomentCoordinates(a_fullCellMap,a_mono,a_refinedCenterDelta);
    }
    else
    {
      moment = changeMomentCoordinates(m_moments,a_mono,a_refinedCenterDelta);
    }

    return moment;
  }

  template <int dim> Real CutCellMoments<dim>::getBdEBMoment(const IvDim             & a_mono,
                                                             const IFData<dim+1>     & a_IFData,
                                                             const IndexTM<Real,dim> & a_refinedCenterDelta)
  {
    // EB moments are 0 if the refined cell (in dim+1) is all in or all out,
    // otherwise they are in the current cell EB moment maps
    Real EBmoment = LARGEREALVAL;

    if (a_IFData.m_allVerticesOut || a_IFData.m_allVerticesIn)
    {
      EBmoment = 0.0;
    }
    else
    {
      EBmoment = changeMomentCoordinates(m_EBmoments,a_mono,a_refinedCenterDelta);
    }

    return EBmoment;
  }

// Methods for reading moments that check the moment exists
  template <int dim> Real CutCellMoments<dim>::getMoment(const IvDim   & a_mono,
                                                         const EBorVol & a_EBorVol) const
  {
    Real moment = LARGEREALVAL;

    if (a_EBorVol == VolMoment)
    {
      // Find the vol in the map
      typename PthMoment::const_iterator it = m_moments.find(a_mono);

      if (it != m_moments.end())
      {
        moment = it->second;
      }
      else
      {
        amrex::Abort("No volume moment in m_moments");
      }
    }
    else if (a_EBorVol == EBMoment)
    {
      // Find the vol in the EBmap
      typename PthMoment::const_iterator it = m_EBmoments.find(a_mono);

      if (it != m_EBmoments.end())
      {
        moment = it->second;
      }
      else
      {
        amrex::Abort("No volume moment in m_bdMoments");
      }
    }
    else
    {
      amrex::Abort("Must ask for EBMoment or VolMoment");
    }

    return moment;
  }

  template <int dim> Real CutCellMoments<dim>::getVol(const EBorVol& a_EBorVol) const
  {
    Real volume;

    // Volume calculated by  integrating x^0 = 1
    IvDim zeroIv = IvDim::Zero;
    volume = getMoment(zeroIv,a_EBorVol);

    return volume;
  }

  template <int dim> IndexTM<Real,dim> CutCellMoments<dim>::getCentroid(const EBorVol& a_EBorVol) const
  {
    // This call checks the string to be either "VolMoment" or "EBMoment"
    Real volume = getVol(a_EBorVol);

    RvDim centroid;

    // Mono= 0,0,...1,0,0
    for (int idir = 0; idir < dim; ++idir)
    {
      IvDim mono     = BASISV_TM<int,dim>(idir);
      centroid[idir] = getMoment(mono,a_EBorVol);

      if (volume < 0.0)
      {
        amrex::Warning("CutCellMoments::getCentroid: Volume fraction is negative");
      }
      else if (volume == 0.0)
      {
        //amrex::Abort("CutCellMoments::getCentroid: Volume fraction is zero");
      }

      centroid[idir] /= volume;
    }

    return centroid;
  }

  template <int dim> Real CutCellMoments<dim>::getResidual(const int & a_iDegree,
                                                           const int & a_normJ) const
  {
    Real retval;
    if (isCovered() || isRegular())
    {
      retval = 0.0;
    }
    else
    {
      retval =  m_residual[a_iDegree][a_normJ];
    }

    return retval;
  }

  template <int dim> void CutCellMoments<dim>::setResidual(const Real& a_value,
                                                           const int & a_iDegree,
                                                           const int & a_normJ)
  {
    if (isCovered() || isRegular())
    {
      amrex::Abort("Invalid assignment to m_residual in covered or full cell");
    }
    else
    {
      m_residual[a_iDegree][a_normJ] = a_value;
    }
  }

  template <int dim> Vector<Real> CutCellMoments<dim>::sliceResidual(const int & a_iDegree) const
  {
    if (isCovered() || isRegular())
    {
      pout()<<"Dim = "<< dim << endl;
      amrex::Abort("Invalid attemtp to slice of m_residual in a covered or full CCM");
    }
    else
    {
      return m_residual[a_iDegree];
    }
  }

  template <int dim> bool CutCellMoments<dim>::isCovered() const
  {
    return m_IFData.m_allVerticesOut;
  }

  template <int dim> bool CutCellMoments<dim>::isRegular() const
  {
    return (m_IFData.m_allVerticesIn && (!m_bdCCOn));
  }

  template <int dim> void CutCellMoments<dim>::print(ostream& a_out) const
  {
    string padding = "  ";
    for (int i = 0; i < GLOBALDIM - dim; i++)
    {
      padding += "  ";
    }

    a_out << padding << "Dim = " << dim << "\n";
    a_out << padding << "\n";

    for (typename PthMoment::const_iterator it = m_moments.begin(); it != m_moments.end(); ++it)
    {
      a_out << padding << "Integral "
            << it->first
            << " = "
            << it->second
            << "\n";

    }

    if (m_moments.size() > 0)
    {
      a_out << padding << "\n";
    }

    for (typename PthMoment::const_iterator it = m_EBmoments.begin(); it != m_EBmoments.end(); ++it)
    {
      std::ios::fmtflags origFlags = a_out.flags();
      int origWidth = a_out.width();
      int origPrecision = a_out.precision();

      a_out << padding << "EBIntegral "
            << it->first
            << " = "
            << it->second
            << "\n";

    }

    if (m_EBmoments.size() > 0)
    {
      a_out << padding << "\n";
    }

    a_out << padding << "Residuals:" << "\n";
    for (int i = 0; i < m_residual.size(); i++)
    {
      a_out << padding << "  " << i << ":";

      std::ios::fmtflags origFlags = a_out.flags();
      int origWidth = a_out.width();
      int origPrecision = a_out.precision();

      for (int j = 0; j < m_residual[i].size(); j++)
      {

        a_out << " "
              << m_residual[i][j];
      }

      a_out << padding << "\n";
    }
    a_out << padding << "\n";

    a_out << padding << "IFData:" << "\n";
    a_out << padding << "\n";
    a_out << m_IFData;

    a_out << padding << "Number of bdCutCellMoments = "
          << m_bdCutCellMoments.size()
          << "\n";
    a_out << padding << "\n";

    for (typename BdCutCellMoments::const_iterator it = m_bdCutCellMoments.begin(); it != m_bdCutCellMoments.end(); ++it)
    {
      a_out << padding << "BdId = " << it->first << "\n";
      a_out << it->second;
    }
  }

  template <int dim> void CutCellMoments<dim>::dump() const
  {
    print(pout());
  }

// Operators
  template <int dim> void CutCellMoments<dim>::operator=(const CutCellMoments<dim> & a_cutCellMoments)
  {
    // Only copy if the objects are distinct
    if (this != &a_cutCellMoments)
    {
      m_moments          = a_cutCellMoments.m_moments;
      m_EBmoments        = a_cutCellMoments.m_EBmoments;
      m_bdCutCellMoments = a_cutCellMoments.m_bdCutCellMoments;
      m_IFData           = a_cutCellMoments.m_IFData;
      m_bdCCOn           = a_cutCellMoments.m_bdCCOn;
      m_residual         = a_cutCellMoments.m_residual;
      m_numActiveBounds  = a_cutCellMoments.m_numActiveBounds;
      m_badNormal        = a_cutCellMoments.m_badNormal;
    }
  }

  template <int dim> ostream& operator<<(ostream                   & a_out,
                                         const CutCellMoments<dim> & a_cutCellMoments)
  {
    a_cutCellMoments.print(a_out);
    return a_out;
  }

}

#endif
