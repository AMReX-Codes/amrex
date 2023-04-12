#include <AMReX_EB_STL_utils.H>
#include <AMReX_EB_triGeomOps_K.H>
#include <AMReX_IntConv.H>
#include <cstring>

namespace amrex
{

namespace {
    // Does line ab intersect with the triangle?
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    bool line_tri_intersects (Real a[3], Real b[3], STLtools::Triangle const& tri)
    {
        if (amrex::max(a[0],b[0]) < amrex::min(tri.v1.x,tri.v2.x,tri.v3.x) ||
            amrex::min(a[0],b[0]) > amrex::max(tri.v1.x,tri.v2.x,tri.v3.x) ||
            amrex::max(a[1],b[1]) < amrex::min(tri.v1.y,tri.v2.y,tri.v3.y) ||
            amrex::min(a[1],b[1]) > amrex::max(tri.v1.y,tri.v2.y,tri.v3.y) ||
            amrex::max(a[2],b[2]) < amrex::min(tri.v1.z,tri.v2.z,tri.v3.z) ||
            amrex::min(a[2],b[2]) > amrex::max(tri.v1.z,tri.v2.z,tri.v3.z))
        {
            return false;
        }
        else
        {
            Real t1[] = {tri.v1.x, tri.v1.y, tri.v1.z};
            Real t2[] = {tri.v2.x, tri.v2.y, tri.v2.z};
            Real t3[] = {tri.v3.x, tri.v3.y, tri.v3.z};
            return 1-tri_geom_ops::lineseg_tri_intersect(a,b,t1,t2,t3);
        }
    }

    // Sign for 3 points on a plane.  This computes the sign of
    // (p2-p1)x(p3-2).  It is used to determine if a point is inside a
    // triangle in 2d.
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    Real sign (Real x1, Real y1, Real x2, Real y2, Real x3, Real y3)
    {
        Real cp = (x2-x1)*(y3-y2) - (x3-x2)*(y2-y1);
        if (std::abs(cp) < std::numeric_limits<Real>::epsilon()) {
            return 0._rt;
        } else {
            return std::copysign(1.0_rt, cp);
        }
    }

    // Does line (x1,y,z)->(x2,y,z) intersect triangle (v1,v2,v3)?
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
    std::pair<bool,Real> edge_tri_intersects (Real x1, Real x2, Real y, Real z,
                                              XDim3 const& v1, XDim3 const& v2,
                                              XDim3 const& v3, XDim3 const& norm,
                                              Real dlevset)
    {
        if ((dlevset > 0._rt && norm.x > 0._rt) || (dlevset < 0._rt && norm.x < 0._rt))
        { // This triangle has the wrong direction // NOLINT(bugprone-branch-clone)
            return std::make_pair(false,0.0_rt);
        }
        else if (x1 > amrex::max(v1.x,v2.x,v3.x) ||
                 x2 < amrex::min(v1.x,v2.x,v3.x) ||
                 y  > amrex::max(v1.y,v2.y,v3.y) ||
                 y  < amrex::min(v1.y,v2.y,v3.y) ||
                 z  > amrex::max(v1.z,v2.z,v3.z) ||
                 z  < amrex::min(v1.z,v2.z,v3.z))
        {
            return std::make_pair(false,0.0_rt);
        }
        else if (norm.x == 0.0_rt)
        {
            // The line is on the triangle's plane or parallel to the plane.
            return std::make_pair(false,0.0_rt);
        }
        else
        {
            Real x0 = (v1.x+v2.x+v3.x) * (1._rt/3._rt);
            Real y0 = (v1.y+v2.y+v3.y) * (1._rt/3._rt);
            Real z0 = (v1.z+v2.z+v3.z) * (1._rt/3._rt);
            Real x = ((norm.x*x0+norm.y*y0+norm.z*z0) - (norm.y*y+norm.z*z)) / norm.x;
            Real s1 = sign(v1.y, v1.z, v2.y, v2.z, y, z);
            Real s2 = sign(v2.y, v2.z, v3.y, v3.z, y, z);
            Real s3 = sign(v3.y, v3.z, v1.y, v1.z, y, z);
            if (s1 == 0._rt || s2 == 0._rt || s3 == 0._rt || (s1 == s2 && s2 == s3)) {
                if (std::abs(x1-x) < std::numeric_limits<Real>::epsilon()) {
                    return std::make_pair(true,x1);
                } else if (std::abs(x2-x) < std::numeric_limits<Real>::epsilon()) {
                    return std::make_pair(true,x2);
                } else if (x>x1 && x<x2) {
                    return std::make_pair(true,x);
                }
            }
            return std::make_pair(false,0.0_rt);
        }
    }
}

void
STLtools::read_stl_file (std::string const& fname, Real scale, Array<Real,3> const& center,
                         int reverse_normal)
{
    if (ParallelDescriptor::IOProcessor()) {
        char header[6];
        header[5] = '\0';
        {
            std::ifstream is(fname, std::istringstream::in|std::ios::binary);
            if (!is.good()) {
                amrex::Abort("STLtools::read_stl_file: failed to open " + fname);
            }
            is.read(header, 5);
        }
        int is_binary = std::strcmp(header, "solid");
        if (is_binary) {
            read_binary_stl_file(fname, scale, center, reverse_normal);
        } else {
            read_ascii_stl_file(fname, scale, center, reverse_normal);
        }
    }

    prepare();
}

void
STLtools::read_binary_stl_file (std::string const& fname, Real scale,
                                Array<Real,3> const& center, int reverse_normal)
{
    if (ParallelDescriptor::IOProcessor()) {
        if (amrex::Verbose()) {
            Print() << "Reading binary STL file "<< fname << "\n";
        }

        IntDescriptor uint32_descr(sizeof(uint32_t), IntDescriptor::ReverseOrder);
        IntDescriptor uint16_descr(sizeof(uint16_t), IntDescriptor::ReverseOrder);
        RealDescriptor real32_descr(FPC::ieee_float, FPC::reverse_float_order, 4);

        std::ifstream is(fname, std::istringstream::in|std::ios::binary);
        if (!is.good()) {
            amrex::Abort("STLtools::read_binary_stl_file: failed to open " + fname);
        }

        char tmp[81];
        tmp[80] = '\0';
        is.read(tmp, 80); // Header - 80 bytes

        uint32_t numtris; // uint32 - Number of triangles - 4 bytes
        amrex::readIntData<uint32_t,uint32_t>(&numtris, 1, is, uint32_descr);
        m_num_tri = static_cast<int>(numtris);
        m_tri_pts_h.resize(m_num_tri);

        if (amrex::Verbose()) {
            Print() << "    Number of triangles: " << m_num_tri << "\n";
        }

        for (int i=0; i < m_num_tri; ++i) {
            is.read(tmp, 50);  // 50 bytes for each triangle. Vertex 1 starts at 12 bytes.
            Real* p = &(m_tri_pts_h[i].v1.x);
            RealDescriptor::convertToNativeFormat(p, 9, tmp+12, real32_descr);
            for (int j = 0; j < 3; ++j) {
                p[0] = p[0] * scale + center[0];
                p[1] = p[1] * scale + center[1];
                p[2] = p[2] * scale + center[2];
                p += 3;
            }
            if (reverse_normal) {
                std::swap(m_tri_pts_h[i].v1, m_tri_pts_h[i].v2);
            }
        }
    }
}

void
STLtools::read_ascii_stl_file (std::string const& fname, Real scale,
                               Array<Real,3> const& center, int reverse_normal)
{
    if (ParallelDescriptor::IOProcessor()) {
        if (amrex::Verbose()) {
            Print() << "Reading binary STL file "<< fname << "\n";
        }

        std::ifstream is(fname, std::istringstream::in);
        if (!is.good()) {
            amrex::Abort("STLtools::read_ascii_stl_file: failed to open " + fname);
        }

        std::string tmp;
        int nlines=0;

        std::getline(is,tmp); //solid <solidname>
        while(!is.eof())
        {
            std::getline(is,tmp);
            if(tmp.find("endsolid")!=std::string::npos)
            {
                break;
            }
            nlines++;
        }

        constexpr int nlines_per_facet=7; //specific to ASCII STLs
        if(nlines%nlines_per_facet!=0)
        {
            Abort("may be there are blank lines in the STL file\n");
        }

        m_num_tri = nlines / nlines_per_facet;
        m_tri_pts_h.resize(m_num_tri);
        static_assert(sizeof(Triangle) == sizeof(Real)*9, "sizeof(Triangle) is wrong");
        Real* p = &(m_tri_pts_h[0].v1.x);

        if (amrex::Verbose()) {
            Print() << "    Number of triangles: " << m_num_tri << "\n";
        }

        is.seekg(0);
        std::getline(is,tmp); //solid <solidname>

        for(int i=0; i < m_num_tri; ++i)
        {
            std::getline(is,tmp); // facet normal
            std::getline(is,tmp); // outer loop

            Real x, y, z;

            for (int iv = 0; iv < 3; ++iv) { // 3 vertices
                is >> tmp >> x >> y >> z;
                *p++ = x * scale + center[0];
                *p++ = y * scale + center[1];
                *p++ = z * scale + center[2];
            }
            std::getline(is,tmp); // read \n

            std::getline(is,tmp); //end loop
            std::getline(is,tmp); //end facet

            if (reverse_normal) {
                std::swap(m_tri_pts_h[i].v1, m_tri_pts_h[i].v2);
            }
        }
    }
}

void
STLtools::prepare ()
{
    ParallelDescriptor::Bcast(&m_num_tri, 1);
    if (!ParallelDescriptor::IOProcessor()) {
        m_tri_pts_h.resize(m_num_tri);
    }
    ParallelDescriptor::Bcast((char*)(m_tri_pts_h.dataPtr()), m_num_tri*sizeof(Triangle));

    //device vectors
    m_tri_pts_d.resize(m_num_tri);
    m_tri_normals_d.resize(m_num_tri);

    Gpu::copyAsync(Gpu::hostToDevice, m_tri_pts_h.begin(), m_tri_pts_h.end(),
                   m_tri_pts_d.begin());

    Triangle const* tri_pts = m_tri_pts_d.data();
    XDim3* tri_norm = m_tri_normals_d.data();

    // Compute normals in case the STL file does not have valid data for normals
    ParallelFor(m_num_tri, [=] AMREX_GPU_DEVICE (int i) noexcept
    {
        Triangle const& tri = tri_pts[i];
        XDim3 vec1{tri.v2.x-tri.v1.x, tri.v2.y-tri.v1.y, tri.v2.z-tri.v1.z};
        XDim3 vec2{tri.v3.x-tri.v2.x, tri.v3.y-tri.v2.y, tri.v3.z-tri.v2.z};
        XDim3 norm{vec1.y*vec2.z-vec1.z*vec2.y,
                   vec1.z*vec2.x-vec1.x*vec2.z,
                   vec1.x*vec2.y-vec1.y*vec2.x};
        Real tmp = 1._rt / std::sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
        tri_norm[i].x = norm.x * tmp;
        tri_norm[i].y = norm.y * tmp;
        tri_norm[i].z = norm.z * tmp;
    });

    ReduceOps<ReduceOpMin,ReduceOpMin,ReduceOpMin,ReduceOpMax,ReduceOpMax,ReduceOpMax> reduce_op;
    ReduceData<Real,Real,Real,Real,Real,Real> reduce_data(reduce_op);
    using ReduceTuple = typename decltype(reduce_data)::Type;
    reduce_op.eval(m_num_tri, reduce_data,
                   [=] AMREX_GPU_DEVICE (int i) -> ReduceTuple
                   {
                       return {amrex::min(tri_pts[i].v1.x,
                                          tri_pts[i].v2.x,
                                          tri_pts[i].v3.x),
                               amrex::min(tri_pts[i].v1.y,
                                          tri_pts[i].v2.y,
                                          tri_pts[i].v3.y),
                               amrex::min(tri_pts[i].v1.z,
                                          tri_pts[i].v2.z,
                                          tri_pts[i].v3.z),
                               amrex::max(tri_pts[i].v1.x,
                                          tri_pts[i].v2.x,
                                          tri_pts[i].v3.x),
                               amrex::max(tri_pts[i].v1.y,
                                          tri_pts[i].v2.y,
                                          tri_pts[i].v3.y),
                               amrex::max(tri_pts[i].v1.z,
                                          tri_pts[i].v2.z,
                                          tri_pts[i].v3.z)};
                   });
    auto const& hv = reduce_data.value(reduce_op);
    m_ptmin.x = amrex::get<0>(hv);
    m_ptmin.y = amrex::get<1>(hv);
    m_ptmin.z = amrex::get<2>(hv);
    m_ptmax.x = amrex::get<3>(hv);
    m_ptmax.y = amrex::get<4>(hv);
    m_ptmax.z = amrex::get<5>(hv);

    if (amrex::Verbose() > 0) {
        amrex::Print() << "    Min: " << m_ptmin << " Max: " << m_ptmax << std::endl;
    }

    // Choose a reference point by extending the normal vector of the first
    // triangle until it's slightly outside the bounding box.
    XDim3 cent0; // centroid of the first triangle
    int is_ref_positive;
    {
        Triangle const& tri = m_tri_pts_h[0];
        cent0 = XDim3{(tri.v1.x + tri.v2.x + tri.v3.x) / 3._rt,
                      (tri.v1.y + tri.v2.y + tri.v3.y) / 3._rt,
                      (tri.v1.z + tri.v2.z + tri.v3.z) / 3._rt};
        // We are computing the normal ourselves in case the stl file does
        // not have valid data on normal.
        XDim3 vec1{tri.v2.x-tri.v1.x, tri.v2.y-tri.v1.y, tri.v2.z-tri.v1.z};
        XDim3 vec2{tri.v3.x-tri.v2.x, tri.v3.y-tri.v2.y, tri.v3.z-tri.v2.z};
        XDim3 norm{vec1.y*vec2.z-vec1.z*vec2.y,
                   vec1.z*vec2.x-vec1.x*vec2.z,
                   vec1.x*vec2.y-vec1.y*vec2.x};
        Real tmp = 1._rt / std::sqrt(norm.x*norm.x + norm.y*norm.y + norm.z*norm.z);
        norm.x *= tmp;
        norm.y *= tmp;
        norm.z *= tmp;
        // Now we need to find out where the normal vector will intersect
        // with the bounding box defined by m_ptmin and m_ptmax.
        Real Lx, Ly, Lz;
        constexpr Real eps = std::numeric_limits<Real>::epsilon();
        if (norm.x > eps) {
            Lx = (m_ptmax.x-cent0.x) / norm.x;
        } else if (norm.x < -eps) {
            Lx = (m_ptmin.x-cent0.x) / norm.x;
        } else {
            Lx = std::numeric_limits<Real>::max();
        }
        if (norm.y > eps) {
            Ly = (m_ptmax.y-cent0.y) / norm.y;
        } else if (norm.y < -eps) {
            Ly = (m_ptmin.y-cent0.y) / norm.y;
        } else {
            Ly = std::numeric_limits<Real>::max();
        }
        if (norm.z > eps) {
            Lz = (m_ptmax.z-cent0.z) / norm.z;
        } else if (norm.z < -eps) {
            Lz = (m_ptmin.z-cent0.z) / norm.z;
        } else {
            Lz = std::numeric_limits<Real>::max();
        }
        Real Lp = std::min({Lx,Ly,Lz});
        if (norm.x > eps) {
            Lx = (m_ptmin.x-cent0.x) / norm.x;
        } else if (norm.x < -eps) {
            Lx = (m_ptmax.x-cent0.x) / norm.x;
        } else {
            Lx = std::numeric_limits<Real>::lowest();
        }
        if (norm.y > eps) {
            Ly = (m_ptmin.y-cent0.y) / norm.y;
        } else if (norm.y < -eps) {
            Ly = (m_ptmax.y-cent0.y) / norm.y;
        } else {
            Ly = std::numeric_limits<Real>::lowest();
        }
        if (norm.z > eps) {
            Lz = (m_ptmin.z-cent0.z) / norm.z;
        } else if (norm.z < -eps) {
            Lz = (m_ptmax.z-cent0.z) / norm.z;
        } else {
            Lz = std::numeric_limits<Real>::lowest();
        }
        if (std::abs(norm.x) < 1.e-5) {
            norm.x = std::copysign(Real(1.e-5), norm.x);
        }
        if (std::abs(norm.y) < 1.e-5) {
            norm.y = std::copysign(Real(1.e-5), norm.y);
        }
        if (std::abs(norm.z) < 1.e-5) {
            norm.z = std::copysign(Real(1.e-5), norm.z);
        }
        Real Lm = std::max({Lx,Ly,Lz});
        Real Leps = std::max(Lp,-Lm) * Real(0.009);
        if (Lp < -Lm) {
            m_ptref.x = cent0.x + (Lp+Leps) * norm.x;
            m_ptref.y = cent0.y + (Lp+Leps) * norm.y;
            m_ptref.z = cent0.z + (Lp+Leps) * norm.z;
            is_ref_positive = true;
        } else {
            m_ptref.x = cent0.x + (Lm-Leps) * norm.x;
            m_ptref.y = cent0.y + (Lm-Leps) * norm.y;
            m_ptref.z = cent0.z + (Lm-Leps) * norm.z;
            is_ref_positive = false;
        }
    }

    // We now need to figure out if the boundary and the reference is
    // outside or inside the object.
    XDim3 ptref = m_ptref;
    int num_isects = Reduce::Sum<int>(m_num_tri, [=] AMREX_GPU_DEVICE (int i) -> int
        {
            if (i == 0) {
                return 1-is_ref_positive;
            } else {
                Real p1[] = {ptref.x, ptref.y, ptref.z};
                Real p2[] = {cent0.x, cent0.y, cent0.z};
                return static_cast<int>(line_tri_intersects(p1, p2, tri_pts[i]));
            }
        });

    m_boundry_is_outside = num_isects % 2 == 0;
}

void
STLtools::fill (MultiFab& mf, IntVect const& nghost, Geometry const& geom,
                Real outside_value, Real inside_value) const
{
    int num_triangles = m_num_tri;

    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();

    const Triangle* tri_pts = m_tri_pts_d.data();
    XDim3 ptmin = m_ptmin;
    XDim3 ptmax = m_ptmax;
    XDim3 ptref = m_ptref;
    Real reference_value = m_boundry_is_outside ? outside_value :  inside_value;
    Real other_value     = m_boundry_is_outside ?  inside_value : outside_value;

    auto const& ma = mf.arrays();

    ParallelFor(mf, nghost, [=] AMREX_GPU_DEVICE (int box_no, int i, int j, int k) noexcept
    {
        Real coords[3];
        coords[0]=plo[0]+static_cast<Real>(i)*dx[0];
        coords[1]=plo[1]+static_cast<Real>(j)*dx[1];
#if (AMREX_SPACEDIM == 2)
        coords[2]=Real(0.);
#else
        coords[2]=plo[2]+static_cast<Real>(k)*dx[2];
#endif
        int num_intersects=0;
        if (coords[0] >= ptmin.x && coords[0] <= ptmax.x &&
            coords[1] >= ptmin.y && coords[1] <= ptmax.y &&
            coords[2] >= ptmin.z && coords[2] <= ptmax.z)
        {
            Real pr[]={ptref.x, ptref.y, ptref.z};
            for (int tr=0; tr < num_triangles; ++tr) {
                if (line_tri_intersects(pr, coords, tri_pts[tr])) {
                    ++num_intersects;
                }
            }
        }
        ma[box_no](i,j,k) = (num_intersects % 2 == 0) ? reference_value : other_value;
    });
    Gpu::streamSynchronize();
}

int
STLtools::getBoxType (Box const& box, Geometry const& geom, RunOn) const
{
    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();

    XDim3 blo{plo[0] + static_cast<Real>(box.smallEnd(0))*dx[0],
              plo[1] + static_cast<Real>(box.smallEnd(1))*dx[1],
#if (AMREX_SPACEDIM == 2)
              0._rt
#else
              plo[2] + static_cast<Real>(box.smallEnd(2))*dx[2]
#endif
    };

    XDim3 bhi{plo[0] + static_cast<Real>(box.bigEnd(0))*dx[0],
              plo[1] + static_cast<Real>(box.bigEnd(1))*dx[1],
#if (AMREX_SPACEDIM == 2)
              0._rt
#else
              plo[2] + static_cast<Real>(box.bigEnd(2))*dx[2]
#endif
    };

    if (blo.x > m_ptmax.x || blo.y > m_ptmax.y || blo.z > m_ptmax.z ||
        bhi.x < m_ptmin.x || bhi.y < m_ptmin.y || bhi.z < m_ptmin.z)
    {
        return m_boundry_is_outside ? allregular : allcovered;
    }
    else
    {
        int num_triangles = m_num_tri;
        const Triangle* tri_pts = m_tri_pts_d.data();
        XDim3 ptmin = m_ptmin;
        XDim3 ptmax = m_ptmax;
        XDim3 ptref = m_ptref;
        int ref_value = m_boundry_is_outside ? 1 : 0;

        ReduceOps<ReduceOpSum> reduce_op;
        ReduceData<int> reduce_data(reduce_op);
        using ReduceTuple = typename decltype(reduce_data)::Type;
        reduce_op.eval(box, reduce_data,
        [=] AMREX_GPU_DEVICE (int i, int j, int k) -> ReduceTuple
        {
            Real coords[3];
            coords[0]=plo[0]+static_cast<Real>(i)*dx[0];
            coords[1]=plo[1]+static_cast<Real>(j)*dx[1];
#if (AMREX_SPACEDIM == 2)
            amrex::ignore_unused(k);
            coords[2]=Real(0.);
#else
            coords[2]=plo[2]+static_cast<Real>(k)*dx[2];
#endif
            int num_intersects=0;
            if (coords[0] >= ptmin.x && coords[0] <= ptmax.x &&
                coords[1] >= ptmin.y && coords[1] <= ptmax.y &&
                coords[2] >= ptmin.z && coords[2] <= ptmax.z)
            {
                Real pr[]={ptref.x, ptref.y, ptref.z};
                for (int tr=0; tr < num_triangles; ++tr) {
                    if (line_tri_intersects(pr, coords, tri_pts[tr])) {
                        ++num_intersects;
                    }
                }
            }

            return (num_intersects % 2 == 0) ? ref_value : 1-ref_value;
        });
        ReduceTuple hv = reduce_data.value(reduce_op);
        Long nfluid = static_cast<Long>(amrex::get<0>(hv));
        Long npts = box.numPts();
        if (nfluid == 0) {
            return allcovered;
        } else if (nfluid == npts) {
            return allregular;
        } else {
            return mixedcells;
        }
    }
}

void
STLtools::fillFab (BaseFab<Real>& levelset, const Geometry& geom, RunOn, Box const&) const
{
    int num_triangles = m_num_tri;

    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();

    const Triangle* tri_pts = m_tri_pts_d.data();
    XDim3 ptmin = m_ptmin;
    XDim3 ptmax = m_ptmax;
    XDim3 ptref = m_ptref;
    Real reference_value = m_boundry_is_outside ? -1.0_rt :  1.0_rt;
    Real other_value     = m_boundry_is_outside ?  1.0_rt : -1.0_rt;

    auto const& a = levelset.array();
    const Box& bx = levelset.box();
    ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
    {
        Real coords[3];
        coords[0]=plo[0]+static_cast<Real>(i)*dx[0];
        coords[1]=plo[1]+static_cast<Real>(j)*dx[1];
#if (AMREX_SPACEDIM == 2)
        coords[2]=Real(0.);
#else
        coords[2]=plo[2]+static_cast<Real>(k)*dx[2];
#endif
        int num_intersects=0;
        if (coords[0] >= ptmin.x && coords[0] <= ptmax.x &&
            coords[1] >= ptmin.y && coords[1] <= ptmax.y &&
            coords[2] >= ptmin.z && coords[2] <= ptmax.z)
        {
            Real pr[]={ptref.x, ptref.y, ptref.z};
            for (int tr=0; tr < num_triangles; ++tr) {
                if (line_tri_intersects(pr, coords, tri_pts[tr])) {
                    ++num_intersects;
                }
            }
        }
        a(i,j,k) = (num_intersects % 2 == 0) ? reference_value : other_value;
    });
}

void
STLtools::getIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                        Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                        Array4<Real const> const& lst ,Geometry const& geom,
                        RunOn, Box const&) const
{
    int num_triangles = m_num_tri;

    const auto plo = geom.ProbLoArray();
    const auto dx  = geom.CellSizeArray();

    const Triangle* tri_pts = m_tri_pts_d.data();
    const XDim3* tri_norm = m_tri_normals_d.data();

    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Array4<Real> const& inter = inter_arr[idim];
        Array4<EB2::Type_t const> const& type = type_arr[idim];
        const Box bx{inter};
        ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept
        {
            Real r = std::numeric_limits<Real>::quiet_NaN();
            if (type(i,j,k) == EB2::Type::irregular) {
                XDim3 p1{plo[0]+static_cast<Real>(i)*dx[0],
                         plo[1]+static_cast<Real>(j)*dx[1],
#if (AMREX_SPACEDIM == 2)
                         Real(0.)
#else
                         plo[2]+static_cast<Real>(k)*dx[2]
#endif
                };
                if (idim == 0) {
                    Real x2 = plo[0]+static_cast<Real>(i+1)*dx[0];
                    int it;
                    for (it=0; it < num_triangles; ++it) {
                        auto const& tri = tri_pts[it];
                        auto tmp = edge_tri_intersects(p1.x, x2, p1.y, p1.z,
                                                       tri.v1, tri.v2, tri.v3,
                                                       tri_norm[it],
                                                       lst(i+1,j,k)-lst(i,j,k));
                        if (tmp.first) {
                            r = tmp.second;
                            break;
                        }
                    }
                    if (it == num_triangles) {
                        r = (lst(i,j,k) > 0._rt) ? p1.x : x2;
                    }
                } else if (idim == 1) {
                    Real y2 = plo[1]+static_cast<Real>(j+1)*dx[1];
                    int it;
                    for (it=0; it < num_triangles; ++it) {
                        auto const& tri = tri_pts[it];
                        auto const& norm = tri_norm[it];
                        auto tmp = edge_tri_intersects(p1.y, y2, p1.z, p1.x,
                                                       {tri.v1.y, tri.v1.z, tri.v1.x},
                                                       {tri.v2.y, tri.v2.z, tri.v2.x},
                                                       {tri.v3.y, tri.v3.z, tri.v3.x},
                                                       {  norm.y,   norm.z,   norm.x},
                                                       lst(i,j+1,k)-lst(i,j,k));
                        if (tmp.first) {
                            r = tmp.second;
                            break;
                        }
                    }
                    if (it == num_triangles) {
                        r = (lst(i,j,k) > 0._rt) ? p1.y : y2;
                    }
                } else {
                    Real z2 = plo[2]+static_cast<Real>(k+1)*dx[2];
                    int it;
                    for (it=0; it < num_triangles; ++it) {
                        auto const& tri = tri_pts[it];
                        auto const& norm = tri_norm[it];
                        auto tmp = edge_tri_intersects(p1.z, z2, p1.x, p1.y,
                                                       {tri.v1.z, tri.v1.x, tri.v1.y},
                                                       {tri.v2.z, tri.v2.x, tri.v2.y},
                                                       {tri.v3.z, tri.v3.x, tri.v3.y},
                                                       {  norm.z,   norm.x,   norm.y},
                                                       lst(i,j,k+1)-lst(i,j,k));
                        if (tmp.first) {
                            r = tmp.second;
                            break;
                        }
                    }
                    if (it == num_triangles) {
                        r = (lst(i,j,k) > 0._rt) ? p1.z : z2;
                    }
                }
            }
            inter(i,j,k) = r;
        });
    }
}

void
STLtools::updateIntercept (Array<Array4<Real>,AMREX_SPACEDIM> const& inter_arr,
                           Array<Array4<EB2::Type_t const>,AMREX_SPACEDIM> const& type_arr,
                           Array4<Real const> const& lst, Geometry const& geom)
{
    auto const& dx = geom.CellSizeArray();
    auto const& problo = geom.ProbLoArray();
    for (int idim = 0; idim < AMREX_SPACEDIM; ++idim) {
        Array4<Real> const& inter = inter_arr[idim];
        Array4<EB2::Type_t const> const& type = type_arr[idim];
        const Box bx{inter};
        amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE (int i, int j, int k)
        {
            if (type(i,j,k) == EB2::Type::irregular) {
                bool is_nan = amrex::isnan(inter(i,j,k));
                if (idim == 0) {
                    if (lst(i,j,k) == Real(0.0) ||
                        (lst(i,j,k) > Real(0.0) && is_nan))
                    {
                        // interp might still be quiet_nan because lst that
                        // was set to zero has been changed by FillBoundary
                        // at periodic bounadries.
                        inter(i,j,k) = problo[0] + static_cast<Real>(i)*dx[0];
                    }
                    else if (lst(i+1,j,k) == Real(0.0) ||
                             (lst(i+1,j,k) > Real(0.0) && is_nan))
                    {
                        inter(i,j,k) = problo[0] + static_cast<Real>(i+1)*dx[0];
                    }
                } else if (idim == 1) {
                    if (lst(i,j,k) == Real(0.0) ||
                        (lst(i,j,k) > Real(0.0) && is_nan))
                    {
                        inter(i,j,k) = problo[1] + static_cast<Real>(j)*dx[1];
                    }
                    else if (lst(i,j+1,k) == Real(0.0) ||
                             (lst(i,j+1,k) > Real(0.0) && is_nan))
                    {
                        inter(i,j,k) = problo[1] + static_cast<Real>(j+1)*dx[1];
                    }
                } else {
                    if (lst(i,j,k) == Real(0.0) ||
                        (lst(i,j,k) > Real(0.0) && is_nan))
                    {
                        inter(i,j,k) = problo[2] + static_cast<Real>(k)*dx[2];
                        }
                    else if (lst(i,j,k+1) == Real(0.0) ||
                             (lst(i,j,k+1) > Real(0.0) && is_nan))
                    {
                        inter(i,j,k) = problo[2] + static_cast<Real>(k+1)*dx[2];
                    }
                }
            }
        });
    }
}

}
