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

#ifndef _FACEINDEX_H_
#define _FACEINDEX_H_

#include "AMReX_LoHiSide.H"
#include "AMReX_VolIndex.H"
#include "AMReX_BaseIndex.H"

namespace amrex
{
  ///
  /**
     FaceIndex is a very lightweight object used
     to distinguish faces.  It has two VolIndex (s)
     and a direction and a way to access them.
  */
  class FaceIndex : public BaseIndex
  {
  public:
    ///
    /**
       if both cellindexes of the vofs are >= 0,
       then the face is interior.  otherwise the face is
       a boundary face.
    */
    FaceIndex(const VolIndex& a_vof1,
              const VolIndex& a_vof2,
              const int& a_direction);
    ///
    /**
       shorthand constructor.
    */
    FaceIndex(const VolIndex& a_vof1,
              const VolIndex& a_vof2);

    ///
    /**
       if both cellindexes of the vofs are >= 0,
       then the face is interior.  otherwise the face is
       a boundary face.
    */
    void define(const VolIndex& a_vof1,
                const VolIndex& a_vof2,
                const int& a_direction);

    ///
    /**
       Shorthand define. First determine if they are neighbors
       then find the direction to pass to full define
    */
    void define(const VolIndex& a_vof1,
                const VolIndex& a_vof2);
    ///
    /**
     */
    FaceIndex();

    ///
    /**
     */
    FaceIndex(const FaceIndex& a_facein);

    ///
    /**
     */
    void define(const FaceIndex& a_facein);

    ///
    /**
     */
    ~FaceIndex();

    ///
    /**
     */
    FaceIndex& operator= (const FaceIndex& a_facein);

    ///
    /**
     */
    bool operator== (const FaceIndex& a_facein) const;

    bool operator< (const FaceIndex& a_rhs) const;

    ///
    /**
     */
    bool operator!= (const FaceIndex& a_facein) const;

    ///
    /**
     */
    const int& direction() const;

    ///
    /**
     */
    const bool& isDefined() const;

    ///
    /**
     */
    const bool& isBoundary() const;

    ///
    /**
       return cell index of the vof on the a_sd side.
       returns -1 if that side of the face is outside
       the domain (boundary face case).
    */
    const int& cellIndex(const Side::LoHiSide& a_sd) const;

    ///
    /**
       return grid index of the vof on the a_sd side.
       may be inside or ouside domain.
    */
    const IntVect& gridIndex(const Side::LoHiSide& a_sd) const;

    ///
    /**
       manufactures the appropriate volindex.
    */
    VolIndex getVoF(const Side::LoHiSide& a_sd) const;

    ///
    /**
       a_vof needs to be one of the two vofs connected by this face
       Standing on a_vof, face opens toward the other vof
       returns 1 if face opens toward positive axis
       returns -1 if face opens toward negative axis
       returns 0 if a_vof is not one of the member vofs.
    */
    int faceSign(const VolIndex& a_vof) const;

    ///for stencil gymnastics
    void shift(const IntVect& a_iv)
    {
      m_loiv += a_iv;
      m_hiiv += a_iv;
    }

  private:

    ///the intvect of the low vof
    /**
     */
    IntVect m_loiv;

    ///the intvect of the high vof
    /**
     */
    IntVect m_hiiv;

    ///
    /**
       The index of the low vof.
       This = -1 if the face is set to be a boundary
       face on the low side of the domain.
    */
    int m_loIndex;

    ///
    /**
       The index of the high vof.
       This = -1 if the face is set to be a boundary
       face on the high side of the domain.
    */
    int m_hiIndex;

    ///direction of the face.
    int m_direction;

    ///
    /**
       true if one of the define functions has been called.
    */
    bool m_isDefined;

    ///
    /**
       true if the face is constructed
       using the boundary face constructor
    */
    bool m_isBoundary;

  };

}
#endif
