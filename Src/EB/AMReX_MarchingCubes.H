#ifndef AMREX_MARCHING_CUBES_H_
#define AMREX_MARCHING_CUBES_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>

#include <map>
#include <memory>
#include <string>

/**
 * \file AMReX_MarchingCubes.H
 *
 * GPU-enabled marching-cubes routines used to extract EB surfaces.
 */

namespace amrex::MC {

//! Initialize internal lookup tables and device buffers for marching cubes.
void Initialize ();
//! Destroy any static storage allocated by Initialize().
void Finalize ();

//! Vertex buffer storing coordinates and normals for extracted surfaces.
struct Vertex
{
    Gpu::DeviceVector<Real> x, y, z, nx, ny, nz;

    //! Resize the vertex buffers to hold \p n entries.
    void resize (int n);

    //! Raw device pointers used when launching kernels.
    GpuArray<Real*,6> dataPtrs () { return {x.data(), y.data(), z.data(), nx.data(), ny.data(), nz.data()}; }
};

//! Triangle index buffer referencing Vertex entries.
struct Triangle
{
    Gpu::DeviceVector<int> v1, v2, v3; // vertex index

    //! Resize the triangle buffer to hold \p n faces.
    void resize (int n);

    //! Raw device pointers for kernel launches.
    GpuArray<int*,3> dataPtrs () { return {v1.data(), v2.data(), v3.data()}; }
};

//! Minimal container holding marching-cubes output for one FAB.
class MCFab
{
public:
    // We will store more data here in the future
    Triangle m_triangles;
    Vertex m_vertices;
};

//! Run marching cubes on signed-distance field \p sdf_fab.
void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab);

//! Write the collected marching-cubes output to an STL file.
void write_stl (std::string const& filename, std::map<int,std::unique_ptr<MCFab>> const& mc_fabs);

}

#endif
