#if (AMREX_SPACEDIM == 3)

#include <limits>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>

#include "AMReX_algoim_integrals.H"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include "algoim_quad.hpp"

namespace {
struct EBshape
{
    template<typename T>
    T operator() (const blitz::TinyVector<T,3>& x) const
    {
        return (x(0) - cent[0])*(norm[0]) + (x(1) - cent[1])*(norm[1]) + (x(2) - cent[2])*(norm[2]);
    }

    template<typename T>
    blitz::TinyVector<T,3> grad(const blitz::TinyVector<T,3>& x) const
    {
        return blitz::TinyVector<double,3>(norm[0],norm[1],norm[2]);
    }

    void setCent(double centx_in, double centy_in, double centz_in)
    {
        cent[0] = centx_in;
        cent[1] = centy_in;
        cent[2] = centz_in;
    }

    void setNormal(double normx_in, double normy_in, double normz_in)
    {
        norm[0] = normx_in;
        norm[1] = normy_in;
        norm[2] = normz_in;
    }

    double cent[3];
    double norm[3];
};
}

namespace amrex 
{

void compute_integrals(MultiFab& intg)
{
    // std::cout << "Algoim computing volfrac " << std::endl;
    // std::cout << std::fixed << std::setprecision(16);

    // This indicates a volume quadrature -- the integral occurs over phi NEGATIVE
    int dim = -1;

    // This isn't used when dim = -1
    int side = -1;

    // This is the degree of the underlying Gaussian quadrature scheme
    int qo = 4;

    const auto& my_factory = dynamic_cast<EBFArrayBoxFactory const&>(intg.Factory());

    const MultiFab*    vfrac = &(my_factory.getVolFrac());
    const MultiCutFab* bcent = &(my_factory.getBndryCent());
//    const MultiCutFab* ccent = &(my_factory.getCentroid());
    const MultiCutFab* bnorm = &(my_factory.getBndryNormal());
    const auto&        flags =   my_factory.getMultiEBCellFlagFab();

    constexpr int i_S_x     =  1-1;
    constexpr int i_S_y     =  2-1;
    constexpr int i_S_z     =  3-1;
    constexpr int i_S_x2    =  4-1;
    constexpr int i_S_y2    =  5-1;
    constexpr int i_S_z2    =  6-1;
    constexpr int i_S_x_y   =  7-1;
    constexpr int i_S_x_z   =  8-1;
    constexpr int i_S_y_z   =  9-1;
    constexpr int i_S_x2_y  = 10-1;
    constexpr int i_S_x2_z  = 11-1;
    constexpr int i_S_x_y2  = 12-1;
    constexpr int i_S_y2_z  = 13-1;
    constexpr int i_S_x_z2  = 14-1;
    constexpr int i_S_y_z2  = 15-1;
    constexpr int i_S_x2_y2 = 16-1;
    constexpr int i_S_x2_z2 = 17-1;
    constexpr int i_S_y2_z2 = 18-1;

    constexpr Real twelfth = 1./12.;
    constexpr Real   offth = 1./144.;

    int ncomp = intg.nComp();

//    Real vol_diff = 0.0;

    for (MFIter mfi(intg,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox();

       auto& gfab = intg[mfi];

       const auto& flag = flags[mfi];
       auto typ = flag.getType(bx);

       blitz::TinyVector<double,3> xmin = {-0.5,-0.5,-0.5};
       blitz::TinyVector<double,3> xmax = { 0.5, 0.5, 0.5};

       if (typ == FabType::covered) {
          gfab.setVal(0.,bx,0,ncomp);
       } else if (typ == FabType::regular) {
          gfab.setVal(0.,bx,0,ncomp);
          gfab.setVal(twelfth,bx,i_S_x2   ,3); // setting S_x^2, S_y^2, S_z^2 = 1/12
          gfab.setVal(  offth,bx,i_S_x2_y2,3); // setting S_x2_y2, S_x2_z2, S_y2_z2 = 1/144
       } else {

          auto const& garr = intg.array(mfi);
          auto const& vfracarr = vfrac->array(mfi);
          auto const& bcentarr = bcent->array(mfi);
          auto const& bnormarr = bnorm->array(mfi);

          auto const& flagarr = flags.array(mfi); 

          const auto lo = amrex::lbound(bx);
          const auto hi = amrex::ubound(bx);

          for (int k = lo.z; k <= hi.z; ++k)
          for (int j = lo.y; j <= hi.y; ++j)
          for (int i = lo.x; i <= hi.x; ++i)
          {
             const Real volfrac = vfracarr(i,j,k);

             // if (volfrac >= (1.0-1.e-12) || volfrac <= 1.0e-12)
             if (flagarr(i,j,k).isRegular())
             {
                 for (int n = 0; n < ncomp; ++n) {
                     garr(i,j,k,n) = 0.0;
                 }
                 garr(i,j,k,i_S_x2) = twelfth;
                 garr(i,j,k,i_S_y2) = twelfth;
                 garr(i,j,k,i_S_z2) = twelfth;
                 garr(i,j,k,i_S_x2_y2) = offth;
                 garr(i,j,k,i_S_x2_z2) = offth;
                 garr(i,j,k,i_S_y2_z2) = offth;
             }
             else if (flagarr(i,j,k).isCovered())
             {
                 for (int n = 0; n < ncomp; ++n) 
                     garr(i,j,k,n) = 0.0;
             } 
             else 
             {
                 const Real centx = bcentarr(i,j,k,0);
                 const Real centy = bcentarr(i,j,k,1);
                 const Real centz = bcentarr(i,j,k,2);
    
                 const Real normx = bnormarr(i,j,k,0);
                 const Real normy = bnormarr(i,j,k,1);
                 const Real normz = bnormarr(i,j,k,2);

                 EBshape eb_phi;
                 eb_phi.setCent(centx,centy,centz);
                 eb_phi.setNormal(normx,normy,normz);

                 // std::cout << "Using centroid: " << centx << " " << centy << " " << centz << std::endl;
                 // std::cout << "Using normal  : " << normx << " " << normy << " " << normz << std::endl;

                 auto q = Algoim::quadGen<3>(eb_phi, Algoim::BoundingBox<double,3>(xmin, xmax), dim, side, qo);
       
                 Real volume   = q([](const auto& x) { return 1.0; });

#if 0
                 if (std::abs(volume - volfrac) > 1.e-12)
                 { 
                    std::cout << "Volume fractions don't match!" << std::endl;
                    std::cout << "VF " << iv << " " << volfrac << " " << volume << std::endl;
                    std::cout << "  Using bndry cent:" << (*bcent)[mfi](iv,0) << " " 
                                                       << (*bcent)[mfi](iv,1) << " " 
                                                       << (*bcent)[mfi](iv,2) << std::endl;
                    std::cout << "  Using bndry norm:" << (*bnorm)[mfi](iv,0) << " " 
                                                       << (*bnorm)[mfi](iv,1) << " " 
                                                       << (*bnorm)[mfi](iv,2) << std::endl;
                    vol_diff = std::max(vol_diff,std::abs(volfrac - volume));
                    // exit(0);
                 } 
#endif
 
                 Real val_S_x   = q([](const auto& x) { return x(0); });
                 Real val_S_y   = q([](const auto& x) { return x(1); });
                 Real val_S_z   = q([](const auto& x) { return x(2); });
                 Real val_S_x2  = q([](const auto& x) { return x(0)*x(0); });
                 Real val_S_y2  = q([](const auto& x) { return x(1)*x(1); });
                 Real val_S_z2  = q([](const auto& x) { return x(2)*x(2); });
 
                 Real val_S_x_y = q([](const auto& x) { return x(0)*x(1); });
                 Real val_S_x_z = q([](const auto& x) { return x(0)*x(2); });
                 Real val_S_y_z = q([](const auto& x) { return x(1)*x(2); });
 
                 Real val_S_x2_y = q([](const auto& x) { return x(0)*x(0)*x(1); });
                 Real val_S_x2_z = q([](const auto& x) { return x(0)*x(0)*x(2); });

                 Real val_S_x_y2 = q([](const auto& x) { return x(0)*x(1)*x(1); });
                 Real val_S_y2_z = q([](const auto& x) { return x(1)*x(1)*x(2); });
 
                 Real val_S_x_z2 = q([](const auto& x) { return x(0)*x(2)*x(2); });
                 Real val_S_y_z2 = q([](const auto& x) { return x(1)*x(2)*x(2); });
 
                 Real val_S_x2_y2 = q([](const auto& x) { return x(0)*x(0)*x(1)*x(1); });
                 Real val_S_x2_z2 = q([](const auto& x) { return x(0)*x(0)*x(2)*x(2); });
                 Real val_S_y2_z2 = q([](const auto& x) { return x(1)*x(1)*x(2)*x(2); });

                 garr(i,j,k,i_S_x) = val_S_x;
                 garr(i,j,k,i_S_y) = val_S_y;
                 garr(i,j,k,i_S_z) = val_S_z;

                 garr(i,j,k,i_S_x2) = val_S_x2;
                 garr(i,j,k,i_S_y2) = val_S_y2;
                 garr(i,j,k,i_S_z2) = val_S_z2;

                 garr(i,j,k,i_S_x_y) = val_S_x_y;
                 garr(i,j,k,i_S_y_z) = val_S_y_z;
                 garr(i,j,k,i_S_x_z) = val_S_x_z;

                 garr(i,j,k,i_S_x2_y) = val_S_x2_y;
                 garr(i,j,k,i_S_x2_z) = val_S_x2_z;

                 garr(i,j,k,i_S_x_y2) = val_S_x_y2;
                 garr(i,j,k,i_S_y2_z) = val_S_y2_z;

                 garr(i,j,k,i_S_x_z2) = val_S_x_z2;
                 garr(i,j,k,i_S_y_z2) = val_S_y_z2;

                 garr(i,j,k,i_S_x2_y2) = val_S_x2_y2;
                 garr(i,j,k,i_S_x2_z2) = val_S_x2_z2;
                 garr(i,j,k,i_S_y2_z2) = val_S_y2_z2;

#if 0
                 auto const& ccentarr = ccent->array(mfi);
                 if (std::abs(val_S_x/volume - (*ccent)[mfi](iv,0)) > 1.e-12 || 
                     std::abs(val_S_y/volume - (*ccent)[mfi](iv,1)) > 1.e-12 || 
                     std::abs(val_S_z/volume - (*ccent)[mfi](iv,2)) > 1.e-12 )
                 { 
                    std::cout << "Centroid doesn't match!" << std::endl;
                    std::cout << "IV                " << iv     << std::endl;
                    std::cout << "VF                " << volume << std::endl;
                    std::cout << "Centroid from amrex:  " << (*ccent)[mfi](iv,0) << " " 
                                                          << (*ccent)[mfi](iv,1) << " " 
                                                          << (*ccent)[mfi](iv,2) << std::endl;

                    std::cout << "Centroid from algoim: " << val_S_x/volume << " " << val_S_y/volume << " " << val_S_z/volume << std::endl;
                    std::cout << "Norm from amrex     : " << (*bnorm)[mfi](iv,0) << " " 
                                                          << (*bnorm)[mfi](iv,1) << " " 
                                                          << (*bnorm)[mfi](iv,2) << "\n" << std::endl;
                    exit(0);
                 }
#endif
             }
          }
       }
    }
#ifdef AMREX_DEBUG
    // If we want to print this, we should do mpi reduce first.
//    amrex::Print() << "Maximum discrepancy in volume fraction " << vol_diff << std::endl;
#endif
}
}
#endif
