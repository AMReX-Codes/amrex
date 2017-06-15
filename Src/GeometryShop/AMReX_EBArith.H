
/*
 *       {_       {__       {__{_______              {__      {__
 *      {_ __     {_ {__   {___{__    {__             {__   {__  
 *     {_  {__    {__ {__ { {__{__    {__     {__      {__ {__   
 *    {__   {__   {__  {__  {__{_ {__       {_   {__     {__     
 *   {______ {__  {__   {_  {__{__  {__    {_____ {__  {__ {__   
 *  {__       {__ {__       {__{__    {__  {_         {__   {__  
 * {__         {__{__       {__{__      {__  {____   {__      {__
 *
 */

#ifndef _EBArith_H___
#define _EBArith_H___
#include <vector>
#include "AMReX_VolIndex.H"
#include "AMReX_EBISBox.H"
#include "AMReX_Stencils.H"
#include "AMReX_IntVectSet.H"
#include <sstream>


namespace amrex
{

/// 
/**
   A collection of useful EB functions.   There is no member data and all functions are static.
   I have tried to trim the (harmongous) Chombo version down to the essentials.
   Here I will put all  the stuff that can be done on a single grid.  EBLevelDataOps will
   hold all the stuff that involves iterating over boxes.
*/
  class EBArith
  {
  public:
    EBArith()
      {   }
    ~EBArith()
      {   }

    static int gradIndex(int ivar, int gradDir)
      {
        return ivar*SpaceDim + gradDir;
      }

    ///
    /**
       Return true if there is a unique adjacent face in the
       given direction.  False if either covered or multivalued.
    */
    static
    bool
    getAdjacentFace(FaceIndex           & a_adjacentFace,
                    const FaceIndex     & a_face,
                    const EBISBox       & a_ebisBox,
                    const Box           & a_domain,
                    const int           & a_idir,
                    const Side::LoHiSide& a_side);

    ///
    static string convertInt(int number)
      {
        stringstream ss;//create a stringstream
        ss << number;//add number to the stream
        return ss.str();//return a string with the contents of the stream
      }

    ///
    static string convertReal(Real number)
      {
        stringstream ss;//create a stringstream
        ss << number;//add number to the stream
        return ss.str();//return a string with the contents of the stream
      }

    ///for hyperbolics, heaven help us all.
    static
    void
    computeCoveredFaces(vector<VolIndex>&     a_coveredFace,
                        IntVectSet&           a_coveredSets,
                        IntVectSet&           a_irregIVS,
                        const int&            a_idir,
                        const Side::LoHiSide& a_sd,
                        const EBISBox&        a_ebisBox,
                        const Box&            a_region);

    ///
    static FaceStencil getInterpStencil(const FaceIndex&     a_face,
                                        const IntVectSet&    a_coarseFineIVS,
                                        const EBISBox&       a_ebisBox,
                                        const Box          & a_domainBox);

    ///
    static void getInterpStencil2D(FaceStencil&          a_sten,
                                   const FaceIndex&      a_face,
                                   const IntVectSet&     a_coarseFineIVS,
                                   const EBISBox&        a_ebisBox,
                                   const Box          &  a_domainBox) ;

    ///
    static void getInterpStencil3D(FaceStencil&         a_sten,
                                   const FaceIndex&     a_face,
                                   const IntVectSet&    a_coarseFineIVS,
                                   const EBISBox&       a_ebisBox,
                                   const Box          & a_domainBox);
    ///
    /** Returns all vofs that can be reached from a_vof via a
        monotone path of length <= than radius
    */
    static void
    getAllVoFsInMonotonePath(std::vector<VolIndex>& a_vofList,
                             const VolIndex&   a_vof,
                             const EBISBox&    a_ebisBox,
                             const int&        a_redistRad);

    ///call the one above this.
    static void
    getAllVoFsInMonotonePath(std::vector<VolIndex>& a_vofList,
                             const IntVect&    a_timesMoved,
                             const IntVect&    a_pathSign,
                             const VolIndex&   a_vof,
                             const EBISBox&    a_ebisBox,
                             const int&        a_redistRad);

    /// get neighbor vofs without regard to path
    static void
    getAllVoFsWithinRadius(std::vector<VolIndex>& a_vofList,
                           const VolIndex&   a_vof,
                           const EBISBox&    a_ebisBox,
                           const int&        a_redistRad);


    ///
    /**
       returns the order of the extrapolation.
       the reason for the last argument is that you might not want the stencil
       to leak over in the noExtrap direction even though you have set a_dist to zero.
       This happens in CF interpolation where you have to really worry about stencil width.
       If orderOfPolynomial = 0,  no first or second derivative stencils included
       If orderOfPolynomial = 1,  no second derivative stencils included
    **/
    static int
    getExtrapolationStencil(VoFStencil&     a_stencil,
                            const RealVect& a_dist,
                            const RealVect& a_dx,
                            const VolIndex& a_startVoF,
                            const EBISBox&  a_ebisBox,
                            int a_orderOfPolynomial = 2,
                            int a_noExtrapThisDirection = -1,
                            IntVectSet*    a_cfivsPtr = NULL,
                            int ivar = 0);
    ///
    static int
    getFirstDerivStencil(VoFStencil&      a_sten,
                         const VolIndex&  a_vof,
                         const EBISBox&   a_ebisBox,
                         const int&       a_idir,
                         const Real&      a_dx,
                         IntVectSet*    a_cfivsPtr,
                         int ivar);

    ///
    static int
    getSecondDerivStencil(VoFStencil&      a_sten,
                          const VolIndex&  a_vof,
                          const EBISBox&   a_ebisBox,
                          const int&       a_idir,
                          const Real&      a_dx,
                          IntVectSet*    a_cfivsPtr,
                          int ivar);

    static int
    getMixedDerivStencil(VoFStencil&      a_sten,
                         const VolIndex&  a_vof,
                         const EBISBox&   a_ebisBox,
                         const int&       a_dir1,
                         const int&       a_dir2,
                         const Real&      a_dx1,
                         const Real&      a_dx2,
                         IntVectSet*    a_cfivsPtr,
                         int ivar);

    static bool
    isVoFHere(VolIndex& a_vof2,
              const std::vector<VolIndex>& a_vofsStencil,
              const IntVect& a_cell2);
    static void
    getLeastSquaresGradSten(VoFStencil&     a_stencil,
                            Real&           a_weight,
                            const RealVect& a_normal  ,
                            const RealVect& a_centroid,
                            const VolIndex& a_vof,
                            const EBISBox&  a_ebisBox,
                            const RealVect& a_dx,
                            const Box& a_domain,
                            int a_ivar);

    static void
    getLeastSquaresGradSten(VoFStencil&     a_stencil,
                            Real&           a_weight,
                            const RealVect& a_normal,
                            const RealVect& a_centroid,
                            const IntVect&  a_quadrant,
                            const VolIndex& a_vof,
                            const EBISBox&  a_ebisBox,
                            const RealVect& a_dx,
                            const Box& a_domain,
                            int a_ivar);

    static void
    johanStencil(bool&               a_dropOrder,
                 std::vector<VoFStencil>& a_pointStencils,
                 std::vector<Real>&       a_distanceAlongLine,
                 const VolIndex&     a_vof,
                 const EBISBox&      a_ebisBox,
                 const RealVect&     a_dx,
                 const IntVectSet&   a_cfivs,
                 int a_ivar=0);

    static void
    johanStencil(bool&               a_dropOrder,
                 std::vector<VoFStencil>& a_pointStencils,
                 std::vector<Real>&       a_distanceAlongLine,
                 const RealVect&     a_normal,
                 const RealVect&     a_bndryCentroid,
                 const VolIndex&     a_vof,
                 const EBISBox&      a_ebisBox,
                 const RealVect&     a_dx,
                 const IntVectSet&   a_cfivs,
                 int a_ivar=0);


    static void
    getLeastSquaresGradStenAllQuad(VoFStencil&          a_stencil,
                                   Real&                a_weight,
                                   const RealVect&      a_normal,
                                   const RealVect&      a_centroid,
                                   const VolIndex&      a_vof,
                                   const EBISBox&       a_ebisBox,
                                   const RealVect&      a_dx,
                                   const Box& a_domain,
                                   int                  a_ivar,
                                   bool                 a_doSymmetric=false);

    static void
    dataRayCast(bool&               a_dropOrder,
                std::vector<VoFStencil>& a_pointStencils,
                std::vector<Real>&       a_distanceAlongLine,
                const RealVect&     a_normal,
                const RealVect&     a_bndryCentroid,
                const VolIndex&     a_vof,
                const EBISBox&      a_ebisBox,
                const RealVect&     a_dx,
                const IntVectSet&   a_cfivs,
                int a_ivar,
                int a_numPoints);

    static void
    calculateWeightingMatrix(RealVect           x0,
                             std::vector<RealVect>&  xp,
                             std::vector<RealVect>&  weightMatrix,
                             bool&              detZero);

    static void
    calculateWeightingMatrixRed(RealVect           x00,
                                std::vector<RealVect>&  xpp,
                                IntVect            dimm,
                                std::vector<RealVect>&  weightMatrix,
                                bool&              deadRed);

    static
    bool
    monotonePathVoFToCellVoF(VolIndex& a_vof2,
                             const VolIndex& a_vof1,
                             const IntVect& a_cell2,
                             const EBISBox& a_ebisBox);


    static bool
    isVoFHere(VolIndex& a_vof2, int& a_whichVoF,
              const std::vector<VolIndex>& a_vofsStencil,
              const IntVect& a_cell2);

    static void
    getVoFsDir(bool& a_hasClose, VolIndex& a_closeVoF,
               bool& a_hasFar,   VolIndex& a_farVoF,
               const EBISBox& a_ebisBox,
               const VolIndex& a_vof,
               int a_idir, Side::LoHiSide a_sd,
               IntVectSet*    a_cfivsPtr);

    ///gets the location in real space of a cell center
    static RealVect  getVofLocation(const VolIndex&  a_vof,
                                    const RealVect&  a_dx,
                                    const RealVect&  a_probLo);

    ///I have misspelled this one time too many
    static RealVect  getVoFLocation(const VolIndex&  a_vof,
                                    const RealVect&  a_dx,
                                    const RealVect&  a_probLo)
      {
        return getVofLocation(a_vof, a_dx, a_probLo);
      }

    ///I have misspelled this one time too many
    static RealVect  getVoFLocation(const VolIndex&  a_vof,
                                    const Real&  a_dx,
                                    const RealVect&  a_probLo)
      {
        RealVect rvdx = a_dx*RealVect::Unit;
        return getVofLocation(a_vof, rvdx, a_probLo);
      }


    ///I have misspelled this one time too many
    static RealVect  getVofLocation(const VolIndex&  a_vof,
                                    const Real&  a_dx,
                                    const RealVect&  a_probLo)
      {
        RealVect rvdx = a_dx*RealVect::Unit;
        return getVofLocation(a_vof, rvdx, a_probLo);
      }


    ///gets the location in real space of a cell center
    static RealVect  getIVLocation(const IntVect&   a_iv,
                                   const RealVect&  a_dx,
                                   const RealVect&  a_probLo);

    ///gets the location in real space of a face center
    static RealVect  getFaceLocation(const FaceIndex& a_face,
                                     const RealVect&  a_dx,
                                     const RealVect&  a_probLo);

    ///gets the location in real space of a face center
    static RealVect  getFaceLocation(const FaceIndex& a_face,
                                     const Real&  a_dx,
                                     const RealVect&  a_probLo)
      {
        RealVect rvdx = a_dx*RealVect::Unit;
        return getFaceLocation(a_face, rvdx, a_probLo);
      }


  };
}
#endif
