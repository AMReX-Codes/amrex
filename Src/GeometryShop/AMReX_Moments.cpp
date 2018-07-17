#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include "AMReX_Moments.H"
#include "AMReX_LSquares.H"
#include "AMReX_RealVect.H"

namespace amrex
{
  /*************/
  // edgeMo Methods //
  void edgeMo::define(const RealVect& a_loPt,
                      const RealVect& a_hiPt,
                      const bool& a_intersectLo,
                      const int& a_direction,
                      const bool& a_covered,
                      const bool& a_regular,
                      const bool& a_dontKnow)
  {
    m_Lo = a_loPt;
    m_Hi = a_hiPt;
    m_intersectLo = a_intersectLo;

    m_varOfInt = a_direction;

    m_edgeLength = m_Hi[m_varOfInt] - m_Lo[m_varOfInt];

    RealVect temp = m_Hi;
    temp += m_Lo;
    temp /= 2.0;

    m_edgeCentroid = temp;

    m_covered  = a_covered;
    m_regular  = a_regular;
    m_dontKnow = a_dontKnow;
  }
  /*************/
  RealVect edgeMo::getLo() const
  {
    return m_Lo;
  }
  /*************/
  RealVect edgeMo::getHi() const
  {
    return m_Hi;
  }
  /*************/
  bool edgeMo::getIntersectLo() const
  {
    return m_intersectLo;
  }
  /*************/

  Real edgeMo::getEdgeLength() const
  {
    return m_edgeLength;
  }
  /*************/

  RealVect edgeMo::getEdgeCentroid() const
  {
    return m_edgeCentroid;
  }

  /*************/
  // access the member data
  bool edgeMo::isCovered() const
  {
    return m_covered;
  }

  bool edgeMo::isRegular() const
  {
    return m_regular;
  }

  bool edgeMo::dontKnow() const
  {
    return m_dontKnow;
  }

  int edgeMo::direction() const
  {
    return m_varOfInt;
  }

  // caclulate moments
  Real edgeMo::moment(const IntVect& a_exponent) const
  {
    IntVect exp = a_exponent;

    exp[m_varOfInt] += 1;

    Real plusOne = exp[m_varOfInt];

    Real lo = 1.0;
    Real hi = 1.0;

    for (int idir = 0; idir < SpaceDim; idir++)
      {
        int num = exp[idir];
        Real curLo = m_Lo[idir];
        Real curHi = m_Hi[idir];

        for (int i = 0; i < num; i++)
          {
            lo *= curLo;
            hi *= curHi;
          }
      }

    return (hi - lo) / plusOne; // plusOne should be at least 1
  }

  // faceMo Methods //
  void faceMo::define(const edgeMo a_edges[4],
                      const int& a_faceNormal,
                      const bool& a_covered,
                      const bool& a_regular,
                      const bool& a_dontKnow)
  {
    for (int i = 0; i < 4; i++)
      {
        m_edges[i] = a_edges[i];
      }

    m_faceNormal = a_faceNormal;

    m_covered  = a_covered;
    m_regular  = a_regular;
    m_dontKnow = a_dontKnow;

    makeNormal();
  }

  // area
  void faceMo::setFaceArea(const Real& a_area)
  {
    m_areaFrac = a_area;
  }

  Real faceMo::getFaceArea() const
  {
    return m_areaFrac;
  }

  // centroid
  void faceMo::setFaceCentroid(const RealVect& a_centroid)
  {
    m_centroid = a_centroid;
  }

  RealVect faceMo::getFaceCentroid() const
  {
    return m_centroid;
  }

  // facetype
  bool faceMo::isRegular() const
  {
    return m_regular;
  }

  bool faceMo::isCovered() const
  {
    return m_covered;
  }

  bool faceMo::dontKnow() const
  {
    return m_dontKnow;
  }

  // bndLength
  Real faceMo::getBdLength() const
  {
    return m_bdLength;
  }

  // normal
  void faceMo::getNormal(Real a_normal[2]) const
  {
    for (int i = 0; i < 2; i++)
      {
        a_normal[i] = m_normalVec[i];
      }
  }

  // edges
  void faceMo::getEdges(edgeMo a_Edges[4]) const
  {
    for (int i = 0; i < 4; i++)
      {
        a_Edges[i] = m_edges[i];
      }
  }

  edgeMo faceMo::retrieveEdge(int& a_iEdge) const
  {
    return m_edges[a_iEdge];
  }

  // faceNormal
  int faceMo::getFaceNormal()const
  {
    return m_faceNormal;
  }

  // construct normal vector
  void faceMo::makeNormal()
  {
    IntVect zeros = IntVect::TheZeroVector();

    for (int index = 0; index < 2; ++index)
      {
        edgeMo HiEdge = m_edges[2*index + 1];
        edgeMo LoEdge = m_edges[2*index];

        m_normalVec[index] = (HiEdge.moment(zeros) - LoEdge.moment(zeros));
      }

    // this sets m_bdlength, too
    normalize(m_normalVec);
  }

  // also set m_bdLength
  void faceMo::normalize(Real a_normalVec[2])
  {
    m_bdLength = 0.0;
    for (int index = 0; index < 2; ++index)
      {
        m_bdLength += a_normalVec[index] * a_normalVec[index];
      }

    if (m_bdLength == 0.0)
      {
        return;
      }

    m_bdLength = sqrt(m_bdLength);
    for (int idir = 0; idir < 2; ++idir)
      {
        a_normalVec[idir] /= m_bdLength;
      }
  }

  // vofMo Methods
  void vofMo::define(const faceMo a_faces[6])
  {
    for (int i = 0; i < 6; ++i)
      {
        edgeMo edges[4];
        a_faces[i].getEdges(edges);

        bool covered  = a_faces[i].isCovered();
        bool regular  = a_faces[i].isRegular();
        bool dontKnow = a_faces[i].dontKnow();

        int faceNormal = a_faces[i].getFaceNormal();

        m_faces[i].define(edges,faceNormal,covered,regular,dontKnow);
      }
    makeNormal();
  }

  void vofMo::makeNormal()
  {
    int order = 0;
    for (int index = 0; index < 3; ++index)
      {
        Moments geom;

        faceMo& HiFace = m_faces[2*index+1];
        faceMo& LoFace = m_faces[2*index];

        Real hiArea = geom.momentCalc2D(order,HiFace)[0];
        Real loArea = geom.momentCalc2D(order,LoFace)[0];

        m_normalVec[index] = hiArea - loArea;
      }

    // this sets m_bdArea, too
    normalize(m_normalVec);
  }

  void vofMo::setNormal(Real a_normalVec[3])
  {
    for (int i = 0; i < 3; i++)
      {
        m_normalVec[i] = a_normalVec[i];
      }
  }

  void  vofMo::getNormal(Real a_normalVec[3]) const
  {
    for (int idir = 0; idir < 3; idir++)
      {
        a_normalVec[idir] = m_normalVec[idir];
      }
  }

  void vofMo::normalize(Real a_normalVec[3])
  {
    m_bdArea = 0.0;
    for (int idir = 0; idir < SpaceDim; ++idir)
      {
        m_bdArea += a_normalVec[idir] * a_normalVec[idir];
      }

    m_bdArea = sqrt(m_bdArea);

    if (m_bdArea != 0.0)
      {
        for (int idir = 0; idir < SpaceDim; ++idir)
          {
            a_normalVec[idir] /= m_bdArea;
          }
      }

    setNormal(a_normalVec);
  }

  void vofMo::getFaces(faceMo a_faces[6]) const
  {
    for (int i = 0; i < 6; ++i)
      {
        edgeMo edges[4];
        m_faces[i].getEdges(edges);

        int faceNormal = m_faces[i].getFaceNormal();

        bool regular  = m_faces[i].isRegular();
        bool covered  = m_faces[i].isCovered();
        bool dontKnow = m_faces[i].dontKnow();

        a_faces[i].define(edges,faceNormal,covered,regular,dontKnow);
      }
  }

  Real vofMo::getBdArea() const
  {
    return m_bdArea;
  }

  // Moments Methods //
  Moments::Moments()
  {
  }

  Vector<Real> Moments::momentCalc3D(const int& a_order,
                                          vofMo&     a_vof)
  {
    Vector<IntVect> list(0),listPlus(0);
    listOfMoments(a_order,list);       // this function knows SpaceDim
    listOfMoments(a_order+1,listPlus);

    Real** A;
    int numRows = listPlus.size()*SpaceDim;
    int numCols = listPlus.size()+list.size();

    LSquares tools;
    tools.allocArray(numRows,numCols,A);

    Real normalVec[3];
    a_vof.getNormal(normalVec);

    Vector<Real> normal(3);
    for (int idir = 0; idir < 3; ++idir)
      {
        normal[idir] = normalVec[idir];
      }

    makeMatrix(list,listPlus,normal,A);

    Vector<Real> rhs(numRows);
    Vector<Real> x(numCols);
    Vector<Real> answer(numRows);
    Vector<Real> answer2DHi;
    Vector<Real> answer2DLo;

    for (int jdir=0;jdir<3;jdir++)
      {

        faceMo faces[6];
        a_vof.getFaces(faces);

        faceMo& hiFace = faces[2*jdir+1];
        faceMo& loFace = faces[2*jdir];

        if (hiFace.isCovered())
          {
            answer2DHi.resize(listPlus.size());

            for (int count = 0; count < listPlus.size(); count++)
              {
                answer2DHi[count] = 0.0;
              }
          }
        else
          {
            answer2DHi = momentCalc2D(a_order+1,hiFace);
          }

        if (loFace.isCovered())
          {
            answer2DLo.resize(listPlus.size());

            for (int count = 0; count < listPlus.size(); count++)
              {
                answer2DLo[count] = 0.0;
              }
          }
        else
          {
            answer2DLo = momentCalc2D(a_order+1,loFace);
          }

        for (int count=0;count<listPlus.size();count++)
          // numRows = length of rhs
          {
            rhs[count*3 + jdir] += answer2DHi[count] -  answer2DLo[count];
          }
      }

    LSquares whichMethod;
    whichMethod.LeastSquares(A,x,rhs);

    tools.freeArray(numRows,numCols,A);

    return x;
  }

  Vector<Real> Moments::momentCalc2D(const int&    a_order,
                                          const faceMo& a_face)
  {
    // 2 vector because 2D, even if SpaceDim=3
    Vector<Real> normal(2);
    Real normalVec[2];

    a_face.getNormal(normalVec);

    for (int j = 0; j < 2; j++)
      {
        normal[j] = (normalVec[j]);
      }

    int faceNormal = a_face.getFaceNormal();

    Vector<IntVect> list(0),listPlus(0);
    listOfMoments(a_order,list); // this function knows SpaceDim
    listOfMoments(a_order+1,listPlus);

    Real** A;
    int numRows = listPlus.size() * 2; // this is a 2 because faces are 2D. If SpaceDim=3 we iterate
    // over SpaceDim and exclude the case where idir = facenormal because the div theorem wouldn't
    // make sense in this context

    int numCols = listPlus.size()+list.size();
    LSquares tools;

    tools.allocArray(numRows,numCols,A);
    int sizeBig = listPlus.size();

    // make the matrix
    makeMatrix(list,listPlus,normal,A,faceNormal);

    Vector<Real> rhs(numRows);

    // x is the unknown
    Vector<Real> x(numCols);

    // make the rhs
    edgeMo edges[4];
    a_face.getEdges(edges);
    Real moment = 0.0;
    int rowIndex = -1;

    for (int exp = 0; exp < sizeBig; ++exp)
      {
        for (int j = 0; j < SpaceDim; j++)
          {
            edgeMo HiEdge;
            edgeMo LoEdge;

            // these little formulae chose which pair of
            // edges is appropriate. No doubt this could be done better
            if (j != faceNormal)
              {
                int low=-1;

                rowIndex += 1;

                if (faceNormal==0)
                  {
                    low = j*j - j;
                  }
                else if (faceNormal==1)
                  {
                    low = j;
                  }
                else if (faceNormal==2)
                  {
                    low = j*j + j;
                  }

                HiEdge = edges[low+1];
                LoEdge = edges[low];

                IntVect monomial = listPlus[exp];
                moment = HiEdge.moment(monomial) - LoEdge.moment(monomial);

                rhs[rowIndex] = moment;
              }
          }
      }

    LSquares whichMethod;
    whichMethod.LeastSquares(A,x,rhs);

    tools.freeArray(numRows,numCols,A);

    return x;
  }

  void  Moments::listOfMoments(const int&       a_order,
                               Vector<IntVect>& a_exponents)
  {
    AMREX_D_TERM(
           for (int i = 0; i <= a_order; ++i),
           for (int j = 0; j <= a_order; ++j),
           for (int k = 0; k <= a_order; ++k))
      if (AMREX_D_TERM(i,+j,+k) == a_order)
        {
          IntVect newexp(AMREX_D_DECL(i,j,k));
          a_exponents.push_back(newexp);
        }
  }

  void Moments::makeMatrix(const Vector<IntVect>& a_list,
                           const Vector<IntVect>& a_listPlus,
                           const Vector<Real>&    a_normalVec,
                           Real**                 a_A,
                           const int&             a_faceNormal)
  {
    int faceNormal; // default argument is
    // absent in SpaceDim = 2 or for a (fully) 3D matrix
    // fully means not merely a 2D face in SpaceDim = 3
    if (SpaceDim == 2)
      {
        faceNormal = 2;
      }
    else
      {
        faceNormal = a_faceNormal;
      }

    int sizeBig = a_listPlus.size();
    int sizeSmall = a_list.size();
    int offset = sizeSmall;
    int rowIndex = -1;
    int vecIndex; // normalVec might only have two components

    for (int i = 0; i < sizeBig; ++i)
      {
        IntVect F = a_listPlus[i];
        for (int j = 0; j < SpaceDim; ++j)
          {
            if (faceNormal == 999)
              {
                vecIndex = j;
              }
            else
              {
                // this gives vecIndex = 0 or 1
                vecIndex = (j*j*(1-faceNormal) + j*(2*faceNormal-1)) / 2;
              }

            if (j != faceNormal)
              {
                int whichCol;

                rowIndex += 1;
                Real coeff = 0.0;

                if (F[j] > 0)
                  {
                    IntVect DivF = F;
                    DivF[j] -= 1;

                    for (int counter = 0; counter < sizeSmall; ++counter)
                      {
                        if (DivF == a_list[counter])
                          {
                            whichCol = counter;
                          }
                      } // this should work as long as there are no duplicates in a_list
                    coeff = F[j];
                    a_A[rowIndex][whichCol] = coeff;
                  }

                a_A[rowIndex][offset+i] = a_normalVec[vecIndex];
              }
          }
      }
  }
} //namespace amrex
