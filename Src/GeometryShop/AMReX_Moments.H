

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

#ifndef _MOMENTS_H_
#define _MOMENTS_H_

#include <AMReX_Array.H>
#include <AMReX_Vector.H>
#include "AMReX_REAL.H"
#include "AMReX_RealVect.H"

namespace amrex
{
  ///
  /**
     This class stores data about an edge of a control volume
  */
  class edgeMo
  {
  public:
    edgeMo()
    {
    }

    void define (const RealVect& a_loPt,
                 const RealVect& a_hiPt,
                 const bool& a_intersectLo,
                 const int& a_direction,
                 const bool& a_covered,
                 const bool& a_regular,
                 const bool& a_dontKnow);

    ///access member data
    RealVect getLo() const;
    RealVect getHi() const;
    bool getIntersectLo() const;

    Real     getEdgeLength() const;
    RealVect getEdgeCentroid() const;

    bool isCovered() const;
    bool isRegular() const;
    bool dontKnow()  const;

    int direction() const;

    ///integrate monomial
    Real moment(const IntVect& a_exponent) const;

  private:
    RealVect m_Lo;
    RealVect m_Hi;

    bool m_intersectLo;

    Real     m_edgeLength;
    RealVect m_edgeCentroid;

    bool m_covered;
    bool m_regular;
    bool m_dontKnow;

    int m_varOfInt;
  };

  /////////////////////////////////class faceMo///////////////////////////////////
  class faceMo
  {
  public:
    faceMo()
    {
    }

    void define(const edgeMo a_edges[4],
                const int&   a_faceNormal,
                const bool&  a_covered,
                const bool&  a_regular,
                const bool&  a_dontKnow);

    void setFaceArea(const Real& a_area);
    Real getFaceArea() const;

    void     setFaceCentroid(const RealVect& a_centroid);
    RealVect getFaceCentroid() const;

    bool isRegular() const;
    bool isCovered() const;
    bool dontKnow()  const;

    Real   getBdLength() const;
    void   getNormal(Real a_normal[2]) const;
    void   getEdges(edgeMo a_Edges[4]) const;
    edgeMo retrieveEdge(int& a_iEdge) const;

    int  getFaceNormal()const;
    void makeNormal();
    void normalize(Real a_normalVec[2]);

  private:
    edgeMo m_edges[4]; //lexigraphic order:0Lo,0Hi,1Lo,1Hi,2Ho,2Hi

    Real m_normalVec[2];

    bool m_covered;
    bool m_regular;
    bool m_dontKnow;

    RealVect m_centroid;
    Real     m_areaFrac;
    int      m_faceNormal;
    Real     m_bdLength;

    // disallowed for all the usual reasons
    void operator=(const faceMo& a_input) = delete;
    faceMo(const faceMo& a_input) = delete;
  };

  /////////////////////////////////class vofMo////////////////////////////////////
  class vofMo
  {
  public:
    vofMo()
    {
    }

    void define(const faceMo a_faces[6]);

    void makeNormal();
    void setNormal(Real a_normalVec[3]);
    void getNormal(Real a_normalVec[3]) const;
    void normalize(Real a_normalVec[3]);

    void getFaces(faceMo a_faces[6]) const;

    Real getBdArea() const;

  private:
    faceMo m_faces[6];

    Real m_normalVec[3];
    Real m_bdArea;
  };

  //////////////////////classMoments//////////////////////////////////////
  class Moments
  {
  public:
    Moments();

    Array<Real> momentCalc3D(const int& a_order,
                                   vofMo&     a_vof);

    Array<Real> momentCalc2D(const int&    a_order,
                                   const faceMo& a_face);

    void  listOfMoments(const int           & a_order,
                        Array<IntVect>& a_exponents);

    void makeMatrix(const Array<IntVect>& a_list,
                    const Array<IntVect>& a_listPlus,
                    const Array<Real>   & a_normalVec,
                    Real**                      a_A,
                    const int&                  a_faceNormal=999
                    ); //use default faceNormal for 3D matrix or when SpaceDim==2
  };


}
#endif
