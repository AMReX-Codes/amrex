#include<AMReX_EB_STL_utils.H>
#include<AMReX_EB_triGeomOps_K.H>

namespace amrex
{
    //================================================================================
    void STLtools::read_ascii_stl_file(std::string fname)
    {
        std::string tmpline,tmp1,tmp2;
        int nlines=0;

        Vector<char> fileCharPtr;
        ParallelDescriptor::ReadAndBcastFile(fname, fileCharPtr);
        std::string fileCharPtrString(fileCharPtr.dataPtr());
        std::istringstream infile(fileCharPtrString, std::istringstream::in);

        if(amrex::Verbose())
            Print()<<"STL file name:"<<fname<<"\n";

        std::getline(infile,tmpline); //solid <solidname>
        while(!infile.eof())
        {
            std::getline(infile,tmpline);
            if(tmpline.find("endsolid")!=std::string::npos)
            {
                break;
            }
            nlines++;
        }

        if(nlines%m_nlines_per_facet!=0)
        {
            Abort("may be there are blank lines in the STL file\n");
        }

        m_num_tri=nlines/m_nlines_per_facet;

        if(amrex::Verbose())
            Print()<<"number of triangles:"<<m_num_tri<<"\n";

        //host vectors
        m_tri_pts_h.resize(m_num_tri*m_ndata_per_tri);
        m_tri_normals_h.resize(m_num_tri*m_ndata_per_normal);

        infile.seekg(0);
        std::getline(infile,tmpline); //solid <solidname>

        for(int i=0;i<m_num_tri;i++)
        {
            std::getline(infile,tmpline);  //facet normal
            std::istringstream fcnormal(tmpline);
            fcnormal>>tmp1>>tmp2
                >>m_tri_normals_h[i*m_ndata_per_normal+0]
                >>m_tri_normals_h[i*m_ndata_per_normal+1]
                >>m_tri_normals_h[i*m_ndata_per_normal+2];

            std::getline(infile,tmpline); // outer loop

            std::getline(infile,tmpline); //vertex 1
            std::istringstream vertex1(tmpline);
            vertex1>>tmp1
                >>m_tri_pts_h[i*m_ndata_per_tri+0]
                >>m_tri_pts_h[i*m_ndata_per_tri+1]
                >>m_tri_pts_h[i*m_ndata_per_tri+2];

            std::getline(infile,tmpline); //vertex 2
            std::istringstream vertex2(tmpline);
            vertex2>>tmp1
                >>m_tri_pts_h[i*m_ndata_per_tri+3]
                >>m_tri_pts_h[i*m_ndata_per_tri+4]
                >>m_tri_pts_h[i*m_ndata_per_tri+5];

            std::getline(infile,tmpline); //vertex 3
            std::istringstream vertex3(tmpline);
            vertex3>>tmp1 //vertex
                >>m_tri_pts_h[i*m_ndata_per_tri+6]
                >>m_tri_pts_h[i*m_ndata_per_tri+7]
                >>m_tri_pts_h[i*m_ndata_per_tri+8];

            std::getline(infile,tmpline); //end loop
            std::getline(infile,tmpline); //end facet
        }

        //device vectors
        m_tri_pts_d.resize(m_num_tri*m_ndata_per_tri);
        m_tri_normals_d.resize(m_num_tri*m_ndata_per_normal);

        Gpu::copy(Gpu::hostToDevice, m_tri_pts_h.begin(),
                m_tri_pts_h.end(), m_tri_pts_d.begin());
        Gpu::copy(Gpu::hostToDevice,
                m_tri_normals_h.begin(), m_tri_normals_h.end(),
                m_tri_normals_d.begin());
    }
    //================================================================================
    void STLtools::stl_to_markerfab(MultiFab& markerfab,Geometry geom,
            Real *point_outside)
    {
        //local variables for lambda capture
        int data_stride   = m_ndata_per_tri;
        int num_triangles = m_num_tri;
        Real outvalue     = m_outside;
        Real invalue      = m_inside;

        const auto plo   = geom.ProbLoArray();
        const auto dx    = geom.CellSizeArray();
        GpuArray<Real,3> outp={point_outside[0],point_outside[1],point_outside[2]};

        const Real *tri_pts=m_tri_pts_d.data();

        for (MFIter mfi(markerfab); mfi.isValid(); ++mfi) // Loop over grids
        {
            const Box& bx = mfi.validbox();
            auto mfab_arr=markerfab[mfi].array();

            ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j, int k)
            {
                Real coords[3],po[3];
                Real t1[3],t2[3],t3[3];

                coords[0]=plo[0]+i*dx[0];
                coords[1]=plo[1]+j*dx[1];
                coords[2]=plo[2]+k*dx[2];

                po[0]=outp[0];
                po[1]=outp[1];
                po[2]=outp[2];

                int num_intersects=0;
                int intersect;

                for(int tr=0;tr<num_triangles;tr++)
                {
                    t1[0]=tri_pts[tr*data_stride+0];
                    t1[1]=tri_pts[tr*data_stride+1];
                    t1[2]=tri_pts[tr*data_stride+2];

                    t2[0]=tri_pts[tr*data_stride+3];
                    t2[1]=tri_pts[tr*data_stride+4];
                    t2[2]=tri_pts[tr*data_stride+5];

                    t3[0]=tri_pts[tr*data_stride+6];
                    t3[1]=tri_pts[tr*data_stride+7];
                    t3[2]=tri_pts[tr*data_stride+8];

                    intersect = tri_geom_ops::lineseg_tri_intersect(po,coords,t1,t2,t3);
                    num_intersects += (1-intersect);
                }
                if(num_intersects%2 == 0)
                {
                    mfab_arr(i,j,k)=outvalue;
                }
                else
                {
                    mfab_arr(i,j,k)=invalue;
                }

            });
        }
    }
    //================================================================================
}
