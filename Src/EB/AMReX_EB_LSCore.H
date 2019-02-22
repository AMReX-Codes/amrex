#ifndef AMReX_EB_LSCore_H_
#define AMReX_EB_LSCore_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AMReX_BCRec.H>
#include <AMReX_Vector.H>
#include <AMReX_AmrCore.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_FluxRegister.H>

#include <AMReX_EB_F.H>
#include <AMReX_EB_utils.H>
#include <AMReX_EB_LSCore_F.H>
#include <AMReX_EB_LSCoreBase.H>
#include <AMReX_EB_levelset.H>

#if defined(BL_USE_SENSEI_INSITU)
namespace amrex {
class AmrMeshInSituBridge;
}
#endif

//using namespace amrex;

namespace amrex {


template<class F>
class LSCore
    : public LSCoreBase
{
public:

    /****************************************************************************
     * public member functions                                                  *
     ***************************************************************************/

    //! Constructor - reads amr in parameters from inputs file
    //!             - sizes multilevel arrays and data structures
    LSCore (EB2::GeometryShop<F> gshop);
    LSCore (EB2::GeometryShop<F> a_gshop,
            const RealBox * rb, int max_level_in, const Vector<int> & n_cell_in, int coord=-1);

    virtual ~LSCore ();

    //! Delete level data. Overrides the pure virtual function in AmrCore
    virtual void ClearLevel (int lev) override;

    //! Make a new level from scratch using provided BoxArray and
    //! DistributionMapping. Only used during initialization. Overrides the pure
    //! virtual function in AmrCore
    virtual void MakeNewLevelFromScratch (int lev, Real time, const BoxArray & ba,
                                          const DistributionMapping & dm) override;

    static const EB2::Level * BuildEBLevel (const EB2::GeometryShop<F> & gshop,
                                            int lev, const Vector<Geometry> & geom );

    static const EB2::Level * BuildEBLevel ( const EB2::GeometryShop<F> gshop,
                                             int lev, const Vector<Geometry> & geom,
                                             int max_crse_level, int required_crse_lev);

private:

    /****************************************************************************
     * private member functions                                                 *
     ***************************************************************************/

    /****************************************************************************
     * private data members                                                     *
     ***************************************************************************/

    //! Used when building finer levels
    EB2::GeometryShop<F> m_gshop;

    /****************************************************************************
     * runtime parameters                                                       *
     ***************************************************************************/

#if defined(BL_USE_SENSEI_INSITU)
    amrex::AmrMeshInSituBridge *insitu_bridge;
#endif
};


#include <AMReX_EB_LSCoreI.H>

}

#endif
