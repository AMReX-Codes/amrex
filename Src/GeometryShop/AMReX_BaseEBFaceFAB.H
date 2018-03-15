#ifndef _BASEEBFACEFAB_H_
#define _BASEEBFACEFAB_H_

#include <cmath>
#include <cstdlib>
#include "AMReX_BaseIFFAB.H"
#include "AMReX_IntVectSet.H"
#include "AMReX_BaseFab.H"
#include "AMReX_EBISBox.H"

namespace amrex
{
       
  ///  Array defined at the Faces of an Box in an EBISBox.
  /**
     An BaseEBFaceFAB is a templated holder for face-centered
     data over a region which consists of the intersection
     of a face-centered box and an EBISBox.
     At every uncovered VoF in this intersection,
     the BaseEBFaceFAB contains a specified number of data values.
     At singly valued faces,
     the data is stored internally in a BaseFab.  At multiply-valued
     faces, the data is stored internally in an  INFab.
     A multiply-valued face is defined to be a face which
     abuts at least one multiply-valued Face.
     BaseEBFaceFAB provides indexing by Face.
     This class does not provide a copy constructor or assignment operator
     for all the usual reasons.
  */
  template <class T>
  class BaseEBFaceFAB
  {
  public:
       
    /// for AggStencil
    int numDataTypes() const
      {
        return 2;
      }
       
    /// for AggStencil
    int   dataType(const BaseIndex& a_baseInd) const
      {
        const FaceIndex* vofPtr = dynamic_cast<const FaceIndex *>(&a_baseInd);
        if (vofPtr == NULL) amrex::Error("Trying to index into an EBCellFAB with something that is not a FaceIndex");
        int retval = 0;
        const Vector<FaceIndex>& irrfaces =m_irrFAB.getFaces();
        for (int iface = 0; iface < irrfaces.size(); iface++)
        {
          if (irrfaces[iface] == *vofPtr)
          {
            retval = 1;
            break;
          }
        }
        return retval;
      }
       
    ///for AggStencil
    long  offset(const BaseIndex& a_baseInd, const int a_ivar) const
      {
        const FaceIndex* vofPtr = dynamic_cast<const FaceIndex *>(&a_baseInd);
        if (vofPtr == NULL) amrex::Error("Trying to index into an EBCellFAB with something that is not a FaceIndex");
        bool irreg = false;
        const Vector<FaceIndex>& irrfaces =m_irrFAB.getFaces();
        for (int iface = 0; iface < irrfaces.size(); iface++)
        {
          if (irrfaces[iface] == *vofPtr)
          {
            irreg = true;
            break;
          }
        }
        long retval = 0;
        if (irreg)
        {
          retval = m_irrFAB.offset(*vofPtr, a_ivar);
        }
        else
        {
          retval = m_regFAB.offset(vofPtr->gridIndex(Side::Hi), a_ivar);
        }
        return retval;
      }
       
    ///for AggStencil
    T* dataPtr(int a_dataType, int a_ivar)
      {
        if (a_dataType == 0)
        {
          return m_regFAB.dataPtr(a_ivar);
        }
        else if (a_dataType == 1)
        {
          return m_irrFAB.dataPtr(a_ivar);
        }
        else
        {
          amrex::Error("bogus datatype input into dataPtr");
        }
        return NULL; //just avoiding compiler warnings with this
      }
       
    const T* dataPtr(int a_dataType, int a_ivar) const
      {
        if (a_dataType == 0)
        {
          return m_regFAB.dataPtr(a_ivar);
        }
        else if (a_dataType == 1)
        {
          return m_irrFAB.dataPtr(a_ivar);
        }
        else
        {
          amrex::Error("bogus datatype input into dataPtr");
        }
        return NULL; //just avoiding compiler warnings with this
       
      }
       
    ///
    BaseEBFaceFAB();
       
    ///
    /**
       Box going into this needs to be cell-centered.
       Data will exist over the surrounding nodes of the box.
    */
    BaseEBFaceFAB(const EBISBox& a_ebisBox,
                  const Box& a_region,
                  int a_iDir, int a_nVar);
       
    ///
    /**
       Box going into this needs to be cell-centered.
       Data will exist over the surrounding nodes of the box.
    */
    virtual void
    define(const EBISBox& a_ebisBox,
           const Box& a_region,
           int a_iDir, int a_nVar);
       
    ///
    virtual ~BaseEBFaceFAB();
       
    ///
    void  clear();
       
    ///
    void setVal(const T& value);
       
    ///
    void setVal(int ivar,const T& value);
       
    ///
    bool isDefined() const;
       
    ///
    int nComp() const ;
       
    ///
    const BaseIFFAB<T>&   getMultiValuedFAB() const;
       
    ///
    BaseIFFAB<T>&  getMultiValuedFAB() ;
       
    ///
    const BaseFab<T>& getSingleValuedFAB() const;
       
    ///
    BaseFab<T>&   getSingleValuedFAB() ;
       
    ///
    int  direction() const;
       
    ///
    /**
       Returns the FACE-CENTERED region of the fab
    */
    const Box& getRegion() const ;
       
    ///
    /**
       Returns the CELL-CENTERED region of the fab
    */
    const Box& getCellRegion() const ;
       
    ///
    T& operator() (const FaceIndex& a_facein, int  a_nVarLoc);
       
    ///
    const T&
    operator() (const FaceIndex& a_facein, int  a_nVarLoc) const;
       
    std::size_t nBytes (const Box& bx, int start_comp, int ncomps) const;

    std::size_t copyToMem (const Box& srcbox,
                           int        srccomp,
                           int        numcomp,
                           void*      dst) const;


    std::size_t copyFromMem (const Box&  dstbox,
                             int         dstcomp,
                             int         numcomp,
                             const void* src);

    /// default false but you can make  make this true for your T if you are sure.
    static bool preAllocatable () { return false; }
    static bool isCopyOMPSafe () { return true; }

    ///
    /**
       Box going into this should be CELL CENTERED.
       copy done over surrounding nodes.
    */
    BaseEBFaceFAB<T>& copy(const BaseEBFaceFAB<T>& source,
                           const Box             & a_srcbox,
                           int                     a_srccomp,
                           const Box&              a_dstbox,
                           int                     a_dstcomp,
                           int                     a_numcomp);
       
       
    BaseIFFAB<T>& plus(const BaseIFFAB<T>& src,
                       const Box&        srcbox,
                       const Box&        destbox,
                       int               srccomp,
                       int               destcomp,
                       int               numcomp=1)
      {
        amrex::Error("it was a bad idea to put this into the general interface");
        return *this;
      }


    /// Invalid but necessary for LevelData<T> to compile
    BaseEBFaceFAB(const Box& a_region, int a_nVar)
      {
        amrex::Error("invalid constructor called for BaseEBFaceFAB");
      }
       
    ///
    const EBISBox& getEBISBox() const;
       
    ///
    void
    setCoveredFaceVal(const T&    a_val,
                      const int&  a_comp,
                      const bool& a_doMulti=true);
       
  protected:
    ///data at faces which abut multi-valued cells
    BaseIFFAB<T> m_irrFAB;
    ///data at faces between two single-valued cells
    BaseFab<T>   m_regFAB;
       
    ///number of data values at each face in BaseEBFaceFAB
    int m_nComp;
       
    //direction of the faces in the BaseEBFaceFAB
    int m_iDir;
       
    ///the face-centered region over which the BaseEBFaceFAB lives
    Box m_regionFace;
       
    ///the cell-centered region over which the BaseEBFaceFAB lives
    Box m_region;
       
    ///has the full define function been called for the BaseEBFaceFAB?
    bool m_isDefined;
       
    EBISBox m_ebisBox;
  private:
    void
    setDefaultValues();

    void operator= (const BaseEBFaceFAB<T>& ebcin);
    BaseEBFaceFAB (const BaseEBFaceFAB<T>& ebcin);

  };
}       

#include "AMReX_BaseEBFaceFABI.H"
       
#endif
