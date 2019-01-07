#if (AMREX_SPACEDIM == 3)

#include <limits>

#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_MultiFabUtil.H>
#include <AMReX_EBMultiFabUtil.H>
#include <AMReX_EB2.H>

#include <AMReX_MLNodeLap_F.H>

#include "AMReX_algoim_integrals.H"

#ifdef _OPENMP
#include <omp.h>
#endif

#include <fstream>
#include "algoim_quad.hpp"

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

namespace amrex 
{

void compute_integrals(MultiFab* intg)
{
    // std::cout << "Algoim computing volfrac " << std::endl;
    // std::cout << std::fixed << std::setprecision(16);

    // This indicates a volume quadrature -- the integral occurs over phi NEGATIVE
    int dim = -1;

    // This isn't used when dim = -1
    int side = -1;

    // This is the degree of the underlying Gaussian quadrature scheme
    int qo = 4;

    EBshape eb_phi;

    const auto& my_factory = dynamic_cast<EBFArrayBoxFactory const&>(intg->Factory());

    const MultiFab*    vfrac = &(my_factory.getVolFrac());
    const MultiCutFab* bcent = &(my_factory.getBndryCent());
    const MultiCutFab* ccent = &(my_factory.getCentroid());
    const MultiCutFab* bnorm = &(my_factory.getBndryNormal());
    const auto&        flags =   my_factory.getMultiEBCellFlagFab();

    int i_S_x     =  1-1;
    int i_S_y     =  2-1;
    int i_S_z     =  3-1;
    int i_S_x2    =  4-1;
    int i_S_y2    =  5-1;
    int i_S_z2    =  6-1;
    int i_S_x_y   =  7-1;
    int i_S_x_z   =  8-1;
    int i_S_y_z   =  9-1;
    int i_S_x2_y  = 10-1;
    int i_S_x2_z  = 11-1;
    int i_S_x_y2  = 12-1;
    int i_S_y2_z  = 13-1;
    int i_S_x_z2  = 14-1;
    int i_S_y_z2  = 15-1;
    int i_S_x2_y2 = 16-1;
    int i_S_x2_z2 = 17-1;
    int i_S_y2_z2 = 18-1;

    // Initialize to the values they would have in a regular cell
    // setVal (val, comp, ncomp, nghost)
    Real twelfth = 1./12.;
    Real   offth = 1./144.;
    intg->setVal(0.     ,i_S_x    ,3,1); // setting S_x, S_y, S_z = 0
    intg->setVal(twelfth,i_S_x2   ,3,1); // setting S_x^2, S_y^2, S_z^2 = 1/12
    intg->setVal(0.     ,i_S_x_y  ,9,1); // setting S_x_y, S_x_z, S_y_z, S_x2_y, S_x2_z, S_x_y2, S_y2_z, S_x_z2, S_y_z2 = 0
    intg->setVal(offth  ,i_S_x2_y2,3,1); // setting S_x2_y2, S_x2_z2, S_y2_z2 = 1/144

    int ncomp = intg->nComp();

    Real vol_diff = 0.0;

    for (MFIter mfi(*intg,true); mfi.isValid(); ++mfi)
    {
       const Box& bx = mfi.growntilebox();

       const int* lo = bx.loVect();
       const int* hi = bx.hiVect();

       auto& gfab = (*intg)[mfi];

       const auto& flag = flags[mfi];
       auto typ = flag.getType(bx);

       int n_count = 0;

       blitz::TinyVector<double,3> xmin = {-0.5,-0.5,-0.5};
       blitz::TinyVector<double,3> xmax = { 0.5, 0.5, 0.5};

       if (typ == FabType::covered) {
           gfab.setVal(0.,bx,0,ncomp);
       } else if (typ == FabType::regular) {
         // Default values have been set above
       } else {

          for (int k = lo[2]; k <= hi[2]; ++k)
          for (int j = lo[1]; j <= hi[1]; ++j)
          for (int i = lo[0]; i <= hi[0]; ++i)
          {
             IntVect iv(i,j,k);

             const Real volfrac = (*vfrac)[mfi](iv);

             if (volfrac < (1.0-1.e-12) && volfrac > 1.0e-12) 
             {
                 const Real centx = (*bcent)[mfi](iv,0);
                 const Real centy = (*bcent)[mfi](iv,1);
                 const Real centz = (*bcent)[mfi](iv,2);
    
                 const Real normx = (*bnorm)[mfi](iv,0);
                 const Real normy = (*bnorm)[mfi](iv,1);
                 const Real normz = (*bnorm)[mfi](iv,2);

                 eb_phi.setCent(centx,centy,centz);
                 eb_phi.setNormal(normx,normy,normz);

                 // std::cout << "Using centroid: " << centx << " " << centy << " " << centz << std::endl;
                 // std::cout << "Using normal  : " << normx << " " << normy << " " << normz << std::endl;

                 auto q = Algoim::quadGen<3>(eb_phi, Algoim::BoundingBox<double,3>(xmin, xmax), dim, side, qo);
       
                 Real volume   = q([](const auto& x) { return 1.0; });

                 if (std::abs(volume - volfrac) > 1.e-12)
                 { 
#if 0
                    std::cout << "Volume fractions don't match!" << std::endl;
                    std::cout << "VF " << iv << " " << volfrac << " " << volume << std::endl;
                    std::cout << "  Using bndry cent:" << (*bcent)[mfi](iv,0) << " " 
                                                       << (*bcent)[mfi](iv,1) << " " 
                                                       << (*bcent)[mfi](iv,2) << std::endl;
                    std::cout << "  Using bndry norm:" << (*bnorm)[mfi](iv,0) << " " 
                                                       << (*bnorm)[mfi](iv,1) << " " 
                                                       << (*bnorm)[mfi](iv,2) << std::endl;
#endif
                    vol_diff = std::max(vol_diff,std::abs(volfrac - volume));
                    // exit(0);
                 } 
 
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

                 gfab(iv,i_S_x) = val_S_x;
                 gfab(iv,i_S_y) = val_S_y;
                 gfab(iv,i_S_z) = val_S_z;

                 gfab(iv,i_S_x2) = val_S_x2;
                 gfab(iv,i_S_y2) = val_S_y2;
                 gfab(iv,i_S_z2) = val_S_z2;

                 gfab(iv,i_S_x_y) = val_S_x_y;
                 gfab(iv,i_S_y_z) = val_S_y_z;
                 gfab(iv,i_S_x_z) = val_S_x_z;

                 gfab(iv,i_S_x2_y) = val_S_x2_y;
                 gfab(iv,i_S_x2_z) = val_S_x2_z;

                 gfab(iv,i_S_x_y2) = val_S_x_y2;
                 gfab(iv,i_S_y2_z) = val_S_y2_z;

                 gfab(iv,i_S_x_z2) = val_S_x_z2;
                 gfab(iv,i_S_y_z2) = val_S_y_z2;

                 gfab(iv,i_S_x2_y2) = val_S_x2_y2;
                 gfab(iv,i_S_x2_z2) = val_S_x2_z2;
                 gfab(iv,i_S_y2_z2) = val_S_y2_z2;

                 n_count++;

#if 0
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
       // std::cout << "Integrated over " << n_count << " cells " << std::endl;
    }
#ifdef AMREX_DEBUG
    amrex::Print() << "Maximum discrepancy in volume fraction " << vol_diff << std::endl;
#endif
}
}
#endif
