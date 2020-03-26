#ifndef WARPX_AVERAGE_H_
#define WARPX_AVERAGE_H_

#include "WarpX.H"

namespace Average{

    using namespace amrex;

    /**
     * \brief Returns the cell-centered average of the floating point data contained
     *        in the input
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *        \c mf_in_arr.
     *
     * \param[in] mf_in_arr floating point data to be averaged
     * \param[in] stag      staggering (index type) of the data
     * \param[in] i         index along x to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] j         index along y to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] k         index along z to access the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \param[in] comp      index along the fourth component of the
     *                      <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                      \c mf_in_arr containing the data to be averaged
     * \return averaged field at cell (i,j,k)
     */
    AMREX_GPU_HOST_DEVICE
    AMREX_FORCE_INLINE
    Real ToCellCenter ( Array4<Real const> const& mf_in_arr,
                        const IntVect stag,
                        const int i,
                        const int j,
                        const int k,
                        const int comp=0 )
    {
        const int sx = stag[0];
        const int sy = stag[1];
#if   (AMREX_SPACEDIM == 2)
        constexpr int sz = 0;
#elif (AMREX_SPACEDIM == 3)
        const int sz = stag[2];
#endif
        return 0.125_rt * (   mf_in_arr(i   ,j   ,k   ,comp)
                            + mf_in_arr(i+sx,j   ,k   ,comp)
                            + mf_in_arr(i   ,j+sy,k   ,comp)
                            + mf_in_arr(i   ,j   ,k+sz,comp)
                            + mf_in_arr(i+sx,j+sy,k   ,comp)
                            + mf_in_arr(i   ,j+sy,k+sz,comp)
                            + mf_in_arr(i+sx,j   ,k+sz,comp)
                            + mf_in_arr(i+sx,j+sy,k+sz,comp) );
    };

    /**
     * \brief Stores the cell-centered average of the floating point data contained
     *        in the input
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *        \c mf_in into the output
     *        <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *        \c mf_out.
     *
     * \param[in,out] mf_out <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       containing the floating point data to be filled with cell-centered averages
     * \param[in]     mf_in  <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       containing the floating point data to be averaged
     * \param[in]     dcomp  offset for the fourth component of the
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       object, extracted from its
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       , where the cell-centered averages will be stored
     * \param[in]     scomp  optional offset for the fourth component of the
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       object, extracted from its
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     *                       , containing the data to be averaged
     * \param[in]     ncomp  optional number of components to loop over for both the
     *                       input and the output
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/structamrex_1_1Array4.html">Array4</a>
     *                       objects, extracted from their respective
     *                       <a href="https://amrex-codes.github.io/amrex/doxygen/classamrex_1_1MultiFab.html">MultiFab</a>
     */
    void ToCellCenter ( MultiFab& mf_out,
                        const MultiFab& mf_in,
                        const int dcomp,
                        const int ngrow,
                        const int scomp=0,
                        const int ncomp=1 );
}

#endif // WARPX_AVERAGE_H_
