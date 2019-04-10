
#ifndef BL_IMULTIFAB_H
#define BL_IMULTIFAB_H

#include <AMReX_BLassert.H>
#include <AMReX_IArrayBox.H>
#include <AMReX_FabArray.H>
#include <AMReX_FabArrayUtility.H>
#include <AMReX_Geometry.H>
#include <memory>

namespace amrex {

/*
* A Collection of IArrayBoxes
*
* The iMultiFab class is publically derived from the
* FabArray<int,IArrayBox> class.  It is a collection (stored as an array) of
* IArrayBoxes useful for storing integer data on a domain defined by
* a union of rectangular regions embedded in a uniform index space.  The
* iMultiFab class extends the function of the underlying FabArray class just
* as the IArrayBox class extends the funtion of BaseFab<int>.  Additional
* member functions are defined for I/O and simple arithmetic operations on
* these aggregate objects.
*
8 This class does NOT provide a copy constructor or assignment operator.
*/
class iMultiFab
    :
    public FabArray<IArrayBox>
{
public:

    /**
    * \brief Constructs an empty iMultiFab.  Data can be defined at a later
    * time using the define member functions inherited
    * from FabArray.
    */
    iMultiFab ();
    /**
    * \brief Constructs a iMultiFab with a valid region defined by bxs and
    * a region of definition defined by the grow factor ngrow.
    */
    iMultiFab (const BoxArray&            bs,
	       const DistributionMapping& dm,
	       int                        ncomp,
	       int                        ngrow,
#ifdef AMREX_STRICT_MODE
               const MFInfo&              info,
               const FabFactory<IArrayBox>& factory);
#else
               const MFInfo&              info = MFInfo(),
               const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>());
#endif

    iMultiFab (const BoxArray&            bs,
	       const DistributionMapping& dm,
	       int                        ncomp,
	       const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
               const MFInfo&              info,
               const FabFactory<IArrayBox>& factory);
#else
               const MFInfo&              info = MFInfo(),
               const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>());
#endif

    /**
    * \brief Make an alias iMultiFab. maketype must be
    * amrex::make_alias.  scomp is the starting component of the
    * alias and ncomp is the number of components in the new aliasing
    * iMultiFab.
    *
    * \param rhs
    * \param maketype
    * \param scomp
    * \param ncomp
    */
    iMultiFab (const iMultiFab& rhs, MakeType maketype, int scomp, int ncomp);

    virtual ~iMultiFab () override = default;

    iMultiFab (iMultiFab&& rhs) noexcept = default;
    iMultiFab& operator= (iMultiFab&& rhs) noexcept = default;

    iMultiFab (const iMultiFab& rhs) = delete;
    iMultiFab& operator= (const iMultiFab& rhs) = delete;

    void operator= (int r);

    /**
    * \brief Returns the minimum value contained in component comp of the
    * iMultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the minimum.  The default is to
    * search only the valid regions of the IArrayBoxes.
    *
    * \param comp
    * \param nghost
    * \param local
    */
    int min (int comp,
             int nghost = 0,
	     bool local = false) const;

    /**
    * \brief Identical to the previous min() function, but confines its
    * search to intersection of Box b and the iMultiFab.
    *
    * \param b
    * \param comp
    * \param nghost
    * \param local
    */
    int min (const Box& b,
             int        comp,
             int        nghost = 0,
	     bool       local = false) const;

    /**
    * \brief Returns the maximum value contained in component comp of the
    * iMultiFab.  The parameter nghost determines the number of
    * boundary cells to search for the maximum.  The default is to
    * search only the valid regions of the IArrayBoxes.
    *
    * \param comp
    * \param nghost
    * \param local
    */
    int max (int comp,
	     int nghost = 0,
	     bool local = false) const;

    /**
    * \brief Identical to the previous max() function, but confines its
    * search to intersection of Box b and the iMultiFab.
    *
    * \param b
    * \param comp
    * \param nghost
    * \param local
    */
    int max (const Box& b,
             int        comp,
             int        nghost = 0,
	     bool       local = false) const;

    /**
    * \brief Returns the sum in component comp
    *
    * \param comp
    * \param nghost
    * \param local
    */
    long sum (int comp, int nghost = 0, bool local = false) const;

    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * specified subregion of the iMultiFab.  The subregion consists
    * of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each IArrayBox in the subregion that should
    * be modified.
    *
    * \param val
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void plus (int val,
               int  comp,
               int  num_comp,
               int  nghost = 0);

    /**
    * \brief Identical to the previous version of plus(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    *
    * \param val
    * \param region
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void plus (int       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);

    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the iMultiFab.  The value
    * of nghost specifies the number of cells in the boundary
    * region that should be modified.
    *
    * \param val
    * \param nghost
    */
    void plus (int val,
               int  nghost);

    /**
    * \brief Adds the scalar value val to the value of each cell in the
    * valid region of each component of the iMultiFab, that also
    * intersects the Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each IArrayBox in
    * the subregion that should be modified.
    *
    * \param val
    * \param region
    * \param nghost
    */
    void plus (int       val,
               const Box& region,
               int        nghost);

    /**
    * \brief Scales the value of each cell in the specified subregion of the
    * iMultiFab by the scalar val (a[i] \<- a[i]*val). The subregion
    * consists of the num_comp components starting at component comp.
    * The value of nghost specifies the number of cells in the
    * boundary region of each IArrayBox in the subregion that should
    * be modified.
    *
    * \param val
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void mult (int val,
               int  comp,
               int  num_comp,
               int  nghost = 0);

    /**
    * \brief Identical to the previous version of mult(), with the
    * restriction that the subregion is further constrained to the
    * intersection with Box region.  The value of nghost specifies the
    * number of cells in the boundary region of each IArrayBox in
    * the subregion that should be modified.
    *
    * \param val
    * \param region
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void mult (int       val,
               const Box& region,
               int        comp,
               int        num_comp,
               int        nghost = 0);

    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the iMultiFab by the scalar val (a[i] \<- a[i]*val).
    * The value of nghost specifies the number of cells in the
    * boundary region that should be modified.
    *
    * \param val
    * \param nghost
    */
    void mult (int val,
               int  nghost = 0);

    /**
    * \brief Scales the value of each cell in the valid region of each
    * component of the iMultiFab by the scalar val (a[i] \<- a[i]*val),
    * that also intersects the Box region.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * IArrayBox in the subregion that should be modified.
    *
    * \param val
    * \param region
    * \param nghost
    */
    void mult (int       val,
               const Box& region,
               int        nghost = 0);

    /**
    * \brief Negates the value of each cell in the specified subregion of
    * the iMultiFab.  The subregion consists of the num_comp
    * components starting at component comp.  The value of nghost
    * specifies the number of cells in the boundary region of each
    * IArrayBox in the subregion that should be modified.
    *
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void negate (int comp,
                 int num_comp,
                 int nghost = 0);

    /**
    * \brief Identical to the previous version of negate(), with the
    * restriction that the subregion is further constrained to
    * the intersection with Box region.
    *
    * \param region
    * \param comp
    * \param num_comp
    * \param nghost
    */
    void negate (const Box& region,
                 int        comp,
                 int        num_comp,
                 int        nghost = 0);

    /**
    * \brief Negates the value of each cell in the valid region of
    * the iMultiFab.  The value of nghost specifies the number of
    * cells in the boundary region that should be modified.
    *
    * \param nghost
    */
    void negate (int nghost = 0);

    /**
    * \brief Negates the value of each cell in the valid region of
    * the iMultiFab that also intersects the Box region.  The value
    * of nghost specifies the number of cells in the boundary region
    * that should be modified.
    *
    * \param region
    * \param nghost
    */
    void negate (const Box& region,
                 int        nghost = 0);

    IntVect minIndex (int comp,
                      int nghost = 0) const;

    IntVect maxIndex (int comp,
                      int nghost = 0) const;

    /**
    * \brief This function adds the values of the cells in mf to the corresponding
    * cells of this iMultiFab.  mf is required to have the same BoxArray or
    * "valid region" as this iMultiFab.  The addition is done only to num_comp
    * components, starting with component number strt_comp.  The parameter
    * nghost specifies the number of boundary cells that will be modified.
    * If nghost == 0, only the valid region of each IArrayBox will be
    * modified.
    *
    * \param mf
    * \param strt_comp
    * \param num_comp
    * \param nghost
    */
    void plus (const iMultiFab& mf,
               int             strt_comp,
               int             num_comp,
               int             nghost);

    /**
    * \brief This function subtracts the values of the cells in mf from the
    * corresponding cells of this iMultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this iMultiFab.  The subtraction is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each IArrayBox will be modified.
    *
    * \param mf
    * \param strt_comp
    * \param num_comp
    * \param nghost
    */
    void minus (const iMultiFab& mf,
                int             strt_comp,
                int             num_comp,
                int             nghost);

    /**
    * \brief This function divides the values of the cells in mf from the
    * corresponding cells of this iMultiFab.  mf is required to have the
    * same BoxArray or "valid region" as this iMultiFab.  The division is
    * done only to num_comp components, starting with component number
    * strt_comp.  The parameter nghost specifies the number of boundary
    * cells that will be modified.  If nghost == 0, only the valid region of
    * each IArrayBox will be modified.  Note, nothing is done to protect
    * against divide by zero.
    *
    * \param mf
    * \param strt_comp
    * \param num_comp
    * \param nghost
    */
    void divide (const iMultiFab& mf,
		 int             strt_comp,
		 int             num_comp,
		 int             nghost);

    /**
    * \brief Add src to dst including nghost ghost cells.
    * The two iMultiFabs MUST have the same underlying BoxArray.
    *
    * \param dst
    * \param src
    * \param srccomp
    * \param dstcomp
    * \param numcomp
    * \param nghost
    */
    static void Add (iMultiFab&       dst,
		     const iMultiFab& src,
		     int             srccomp,
		     int             dstcomp,
		     int             numcomp,
		     int             nghost);

    /**
    * \brief Copy from src to dst including nghost ghost cells.
    * The two iMultiFabs MUST have the same underlying BoxArray.
    *
    * \param dst
    * \param src
    * \param srccomp
    * \param dstcomp
    * \param numcomp
    * \param nghost
    */
    static void Copy (iMultiFab&       dst,
                      const iMultiFab& src,
                      int             srccomp,
                      int             dstcomp,
                      int             numcomp,
                      int             nghost);

    /**
    * \brief Subtract src from dst including nghost ghost cells.
    * The two iMultiFabs MUST have the same underlying BoxArray.
    *
    * \param dst
    * \param src
    * \param srccomp
    * \param dstcomp
    * \param numcomp
    * \param nghost
    */
    static void Subtract (iMultiFab&       dst,
			  const iMultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    /**
    * \brief Multiply dst by src including nghost ghost cells.
    * The two iMultiFabs MUST have the same underlying BoxArray.
    *
    * \param dst
    * \param src
    * \param srccomp
    * \param dstcomp
    * \param numcomp
    * \param nghost
    */
    static void Multiply (iMultiFab&       dst,
			  const iMultiFab& src,
			  int             srccomp,
			  int             dstcomp,
			  int             numcomp,
			  int             nghost);

    /**
    * \brief Divide dst by src including nghost ghost cells.
    * The two iMultiFabs MUST have the same underlying BoxArray.
    *
    * \param dst
    * \param src
    * \param srccomp
    * \param dstcomp
    * \param numcomp
    * \param nghost
    */
    static void Divide (iMultiFab&       dst,
			const iMultiFab& src,
			int             srccomp,
			int             dstcomp,
			int             numcomp,
			int             nghost);

    virtual void define (const BoxArray&            bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 const IntVect&             ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<IArrayBox>& factory) override;
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>()) override;
#endif

    virtual void define (const BoxArray&            bxs,
			 const DistributionMapping& dm,
			 int                        nvar,
			 int                        ngrow,
#ifdef AMREX_STRICT_MODE
			 const MFInfo&              info,
                         const FabFactory<IArrayBox>& factory) override;
#else
			 const MFInfo&              info = MFInfo(),
                         const FabFactory<IArrayBox>& factory = DefaultFabFactory<IArrayBox>()) override;
#endif

    static void Initialize ();
    static void Finalize ();
};

std::unique_ptr<iMultiFab>
OwnerMask (FabArrayBase const& mf, const Periodicity& period);

}

#endif /*BL_IMULTIFAB_H*/
