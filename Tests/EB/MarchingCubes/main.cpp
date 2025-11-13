
#include <AMReX.H>
#include <AMReX_EB2.H>
#include <AMReX_ParmParse.H>
#include <AMReX_PlotFileUtil.H>

using namespace amrex;

void main_main ()
{
    int nx = 128;
    int ny = 128;
    int nz = 128;
    int max_grid_size = 64;
    std::string plot_file{"plt"};
    int which_stl = 0;
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
        pp.query("plot_file", plot_file);
        pp.query("which_stl", which_stl);

        ParmParse ppeb2("eb2");
        std::string geom_type("stl");
        ppeb2.add("geom_type", geom_type);
    }

    {
        std::string stl_file;
        Real stl_scale = 1.0;
        std::vector<Real> stl_center{0.0, 0.0, 0.0};

        if (which_stl == 0) {
            stl_file = "AileM6_with_sharp_TE.stl";
            nx = ny = nz = 256;
            stl_scale = 0.001;
            stl_center = std::vector<Real>{-0.8, -0.8, 0.};
        } else if (which_stl == 1) {
            stl_file = "sphere.stl";
        } else if (which_stl == 2) {
            stl_file = "cube.stl";
        } else if (which_stl == 3) {
            stl_file = "Sphericon.stl";
            stl_scale = 0.001;
        } else if (which_stl == 4) {
            stl_file = "stls/blindt.stl";
            stl_scale = 0.012;
            stl_center = std::vector<Real>{-1., 0., -1.};
        } else if (which_stl == 5) {
            stl_file = "stls/cone_2.stl";
            stl_center = std::vector<Real>{-1., -1., -0.45};
        } else if (which_stl == 6) {
            stl_file = "stls/cone.stl";
            stl_center = std::vector<Real>{-1., -1., -0.25};
        } else if (which_stl == 7) {
            stl_file = "stls/cube.stl";
            stl_scale = 2.;
            stl_center = std::vector<Real>{-1., -1., -1.};
        } else if (which_stl == 8) {
            stl_file = "stls/cylinder.stl";
            stl_scale = 0.4;
            stl_center = std::vector<Real>{0., 0., -1.};
        } else if (which_stl == 9) {
            stl_file = "stls/elbow1.stl";
            stl_scale = 1.7;
            stl_center = std::vector<Real>{1., 0., 1.};
        } else if (which_stl == 10) {
            stl_file = "stls/elbow2.stl";
            stl_scale = 1.0;
            stl_center = std::vector<Real>{1., 0., 1.};
        } else if (which_stl == 11) {
            stl_file = "stls/elbow3.stl";
            stl_scale = 1.0;
            stl_center = std::vector<Real>{1., 0., 1.};
        } else if (which_stl == 12) {
            stl_file = "stls/elbow5.stl";
            stl_scale = 0.35;
            stl_center = std::vector<Real>{1., 0., 1.};
        } else if (which_stl == 13) {
            stl_file = "stls/elbow.stl";
            stl_scale = 0.022;
            stl_center = std::vector<Real>{-1., -1., 0.};
        } else if (which_stl == 14) {
            stl_file = "stls/halfsphere.stl";
            stl_scale = 0.16;
            stl_center = std::vector<Real>{0., 0., 0.49};
        } else if (which_stl == 15) {
            // VERY slow
            stl_file = "stls/pipe_area_change1.stl";
            stl_scale = 0.0155;
            stl_center = std::vector<Real>{0., 0., -0.1};
        } else if (which_stl == 16) {
            // VERY slow
            stl_file = "stls/pipe_area_change2.stl";
            stl_scale = 0.012;
            stl_center = std::vector<Real>{0., 0., 0.1};
        } else if (which_stl == 17) {
            stl_file = "stls/pipe_area_change3.stl";
            stl_scale = 0.018;
            stl_center = std::vector<Real>{0., -1., 0.};
        } else if (which_stl == 18) {
            stl_file = "stls/rect.stl";
            stl_scale = 6.6;
            stl_center = std::vector<Real>{-0.5, -0.5, -1.};
        } else if (which_stl == 19) {
            // VERY slow
            stl_file = "stls/unit_sphere.stl";
            stl_scale = 4.;
            stl_center = std::vector<Real>{0., 0., 0.};
        } else if (which_stl == 20) {
            stl_file = "stls/stl-annulus/CA_500.stl";
            stl_scale = 0.002;
            stl_center = std::vector<Real>{0., 0., -1.};
        } else if (which_stl == 21) {
            stl_file = "stls/stl-annulus/CA_600.stl";
            stl_scale = 0.002;
            stl_center = std::vector<Real>{0., 0., -1.};
        } else if (which_stl == 22) {
            stl_file = "stls/stl-annulus/CA_950.stl";
            stl_scale = 0.002;
            stl_center = std::vector<Real>{0., 0., -1.};
        } else if (which_stl == 23) {
            stl_file = "stls/stl-annulus/CA_ALL.stl";
            stl_scale = 0.002;
            stl_center = std::vector<Real>{0., 0., -1.};
        } else if (which_stl == 24) {
            stl_file = "AileM6_with_sharp_TE.stl";
            xmin =  -40.;
            xmax = 1280.;
            ymin =  -40.;
            ymax = 1280.;
            zmin = -40.;
            zmax =  40.;
            nx = 528;
            ny = 528;
            nz =  32;
        } else if (which_stl == 25) {
            stl_file = "stls/HLPW-4_CRM-HL_40-37_Nominal_v2.stl";
            stl_scale = 1.5e-5;
            stl_center = std::vector<Real>{-0.5, 0.01, 0.};
            xmin = -0.5;
            xmax =  0.5;
            ymin =  0.0;
            ymax =  0.5;
            zmin =  0.0;
            zmax =  0.25;
            nx = 512;
            ny = 256;
            nz = 128;
        } else if (which_stl == 26) {
            stl_file = "armadillo.stl";
            xmin = -64.;
            xmax =  64.;
            ymin = -64.;
            ymax = 128.;
            zmin = -64.;
            zmax =  64;
            nx = 128;
            ny = 192;
            nz = 128;
        }

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

#if 0

    auto const& factory = makeEBFabFactory(geom, ba, dm, {1,1,1}, EBSupport::full);
    MultiFab const& vfrc = factory->getVolFrac();
    amrex::VisMF::Write(vfrc, "vfrc");
//#if 0
    amrex::WriteMLMF(plot_file, {&vfrc}, {geom});
    MultiFab foo(vfrc.boxArray(), vfrc.DistributionMap(), 2, vfrc.nGrowVect());
    foo.setVal(1.0);
    MultiFab::Copy(foo, vfrc, 0, 0, 1, vfrc.nGrowVect());
    MultiFab::Subtract(foo, vfrc, 0, 1, 1, vfrc.nGrowVect());
    amrex::WriteMLMF("plt2", {&foo}, {geom});
//#endif
#endif
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    main_main();
    amrex::Finalize();
}
