#ifndef _STENCILS_H_
#define _STENCILS_H_

#include "AMReX_VolIndex.H"
#include "AMReX_BaseIndex.H"
#include "AMReX_FaceIndex.H"
#include "AMReX_REAL.H"


namespace amrex
{
  ///
  /**
     class for vofstencil and basestencil to inherit from
     so AggStencil can unify them.
  */
  class BaseStencil
  {
  public:
    ///
    BaseStencil()
    {
    }

    ///
    virtual ~BaseStencil()
    {
    }

    ///
    virtual int size() const = 0;

    ///
    virtual const BaseIndex& index(int isten) const = 0;


    ///
    virtual const Real& weight(int isten) const = 0;

    ///
    virtual const int& variable(int isten) const = 0;
  };

  /// VoF-centered stencil
  /**
     This stencil is a fundamental tool for building
     eb applications by making the ability to cache
     stencil information. This object consists of
     a vector of VoFs and a vector of weights.
  */
  class VoFStencil: public BaseStencil
  {
  public:
    ///
    /**
       default constructor.  makes empty vectors.
    */
    VoFStencil();

    ///
    virtual ~VoFStencil();

    ///
    void clear();

    ///
    /**
       add a VoF to the Stencil, with it's associated weight
       it is required that all VoFs are in the same EBIndexSpace
       If the vof is already in the
       stencil, add the weights.
    */
    void
    add(const VolIndex& vof,const Real& weight, int ivar = 0);

    ///
    /**
       number of VoFs in the Stencil
    */
    virtual int size() const;

    ///
    /**
       access a VoF
    */
    const VolIndex&  vof(int i) const;

 

    ///
    virtual BaseIndex& index(int i) const
    {
      return (BaseIndex&)(vof(i));
    }

    ///
    /**
       access a weight
    */
    virtual const Real& weight(int i) const;

    ///
    /**
       access a weight
    */
    virtual Real& weight(int i);

    ///
    /**
       shift all entries by a_shift
    */
    void shift(const IntVect& a_shift);

    ///
    virtual const int& variable(int i) const;

    ///
    virtual int& variable(int i);

    ///
    /**
       add all faces and weights of inputs
       to this.  If a vof is already in the
       stencil, add the weights.
       only addition is well-defined here
       as far as arithmatic operations are
       concerned.
    */
    VoFStencil&
    operator+=(const VoFStencil& a_vofsten);

    ///
    /**
     */
    void operator*=(const Real& scaling);

    ///
    VoFStencil&
    operator=(const VoFStencil& a_vofsten);

    ///
    VoFStencil(const VoFStencil& a_vofstenin);

    void setAllVariables(int a_var)
    {
      if (vofs.size() > 0)
        {
          variables = Vector<int>(vofs.size(), a_var);
        }
    }

  void outputToPout() const;

  protected:

    /// the VoFs
    Vector<VolIndex> vofs;
    /// the weights
    Vector<Real> weights;

    Vector<int> variables;
  };

  ///
  /**
     Face-centered Stencil for embedded boundary applications.
  */
  class FaceStencil: public BaseStencil
  {
  public:
    ///
    FaceStencil();

    ///
    virtual ~FaceStencil();

    ///
    void clear();

    ///
    /**
       add an Face and it's weight
       If the face is already in the
       stencil, add the weights.
    */
    void  add(const FaceIndex& face,const Real& weight, int variable=0);

    ///
    /**
       shift all entries by a_shift
    */
    void shift(const IntVect& a_shift);

    ///
    /**
       number of Faces in the Stencil
    */
    virtual int   size() const;

    ///
    void setAllVariables(int a_var)
    {
      if (faces.size() > 0)
        {
          variables = Vector<int>(faces.size(), a_var);
        }
    }

    ///
    /** access an Face
     */
    const FaceIndex&  face(int i) const;

    ///
    /**
       access a weight
    */
    virtual const Real& weight(int i) const;

    ///
    virtual BaseIndex& index(int i) const
    {
      return (BaseIndex&)(face(i));
    }

    ///
    virtual const int& variable(int i) const;

    ///
    virtual int& variable(int i);

    ///
    /**
       add all faces and weights of inputs
       to this.  If a face is already in the
       stencil, add the weights.
       only addition is well-defined here
       as far as arithmatic operations are
       concerned.
    */
    FaceStencil& operator+=(const FaceStencil& a_facesten);

    ///
    FaceStencil& operator=(const FaceStencil& a_facesten);

    ///
    FaceStencil(const FaceStencil& a_facesten);

    ///
    /**
     */
    void operator*=(const Real& scaling);

    void outputToPout() const;
  private:

    /// the Faces
    Vector<FaceIndex> faces;
    /// the weights
    Vector<Real> weights;
    /// the variable numbers
    Vector<int> variables;

  };

  class EBCellFAB;

  ///the variables numbers live in the stencil
  extern Real applyVoFStencil(const VoFStencil& a_vofSten, const EBCellFAB& a_fab);

}




#endif
