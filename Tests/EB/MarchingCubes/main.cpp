
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void main_main ()
{
    int nx = 64;
    int ny = 64;
    int nz = 64;
    int max_grid_size = 32;
    Real xmin = -1.2;
    Real xmax =  1.2;
    Real ymin = -1.2;
    Real ymax =  1.2;
    Real zmin = -1.2;
    Real zmax =  1.2;
    {
        ParmParse pp;
        pp.query("nx", nx);
        pp.query("ny", ny);
        pp.query("nz", nz);
        pp.query("max_grid_size", max_grid_size);

        ParmParse ppeb2("eb2");
        std::string geom_type("stl");
        ppeb2.add("geom_type", geom_type);
        ppeb2.add("test_marching_cubes", 1);
    }

    {
        std::string stl_file("cube.stl");
        Real stl_scale = 1.0;
        std::vector<Real> stl_center{0.0, 0.0, 0.0};

        ParmParse pp("eb2");
        pp.add("stl_file", stl_file);
        pp.add("stl_scale", stl_scale);
        pp.addarr("stl_center", stl_center);
    }

    Geometry geom(Box(IntVect(0),IntVect(AMREX_D_DECL(nx,ny,nz))),
                  RealBox({AMREX_D_DECL(xmin,ymin,zmin)},
                          {AMREX_D_DECL(xmax,ymax,zmax)}),
                  0, {AMREX_D_DECL(0,0,0)});
    BoxArray ba(geom.Domain());
    ba.maxSize(max_grid_size);
    DistributionMapping dm(ba);

    double t0 = amrex::second();
    EB2::BuildMultiValuedMultiCut(geom,0,0);
    double t1 = amrex::second();
    amrex::Print() << "Build time: " << t1-t0 << "\n";
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}
