#ifndef AMREX_EB_STL_UTILS_H_
#define AMREX_EB_STL_UTILS_H_

#include <AMReX_Config.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_MultiFab.H>
#include <AMReX_REAL.H>
#include <AMReX_Vector.H>
#include <AMReX_Geometry.H>
#include <AMReX_FArrayBox.H>
#include <AMReX_Box.H>

namespace amrex
{
    class STLtools
    {
        private:

            //host vectors
            Gpu::PinnedVector<Real> m_tri_pts_h;
            Gpu::PinnedVector<Real> m_tri_normals_h;

            //device vectors
            Gpu::DeviceVector<amrex::Real> m_tri_pts_d;
            Gpu::DeviceVector<amrex::Real> m_tri_normals_d;

            int  m_num_tri=0;
            int  m_ndata_per_tri=9;    //three points x 3 coordinates
            int  m_ndata_per_normal=3; //three components
            int  m_nlines_per_facet=7; //specific to ASCII STLs
            Real m_inside  = -1.0;
            Real m_outside =  1.0;


        public:

            void read_ascii_stl_file(std::string fname);
            void stl_to_markerfab(MultiFab& markerfab,
                    Geometry geom,Real *point_outside);

    };
}
#endif
