#ifndef AMREX_MARCHING_CUBES_H_
#define AMREX_MARCHING_CUBES_H_

#include <AMReX_FArrayBox.H>
#include <AMReX_Geometry.H>
#include <AMReX_GpuContainers.H>

#include <map>
#include <memory>
#include <string>

namespace amrex::MC {

void Initialize ();
void Finalize ();

struct Vertex
{
    Gpu::DeviceVector<Real> x, y, z, nx, ny, nz;

    void resize (int n);

    GpuArray<Real*,6> dataPtrs () { return {x.data(), y.data(), z.data(), nx.data(), ny.data(), nz.data()}; }
};

struct Triangle
{
    Gpu::DeviceVector<int> v1, v2, v3; // vertex index

    void resize (int n);

    GpuArray<int*,3> dataPtrs () { return {v1.data(), v2.data(), v3.data()}; }
};

class MCFab
{
public:
    // We will store more data here in the future
    Triangle m_triangles;
    Vertex m_vertices;
};

void marching_cubes (Geometry const& geom, FArrayBox& sdf_fab, MCFab& mc_fab);

void write_stl (std::string const& filename, std::map<int,std::unique_ptr<MCFab>> const& mc_fabs);

}

#endif
