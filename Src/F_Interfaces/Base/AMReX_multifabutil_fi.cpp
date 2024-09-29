#include <AMReX_MultiFabUtil.H>

using namespace amrex;

extern "C"
{
    void amrex_fi_average_down (const MultiFab* S_fine, MultiFab* S_crse,
                                const Geometry* fgeom, const Geometry* cgeom,
                                int scomp, int ncomp, int rr)
    {
        amrex::average_down(*S_fine, *S_crse, *fgeom, *cgeom, scomp, ncomp, rr);
    }

    void amrex_fi_average_down_cell_node (const MultiFab* S_fine, MultiFab* S_crse,
                                          int scomp, int ncomp, int rr)
    {
        amrex::average_down(*S_fine, *S_crse, scomp, ncomp, rr);
    }

    void amrex_fi_average_down_faces (MultiFab const* fmf[], MultiFab* cmf[],
                                      Geometry const* cgeom, int scomp, int ncomp,
                                      int rr)
    {
        Array<MultiFab,AMREX_SPACEDIM> fine
            {AMREX_D_DECL(MultiFab(*fmf[0], amrex::make_alias, scomp, ncomp),
                          MultiFab(*fmf[1], amrex::make_alias, scomp, ncomp),
                          MultiFab(*fmf[2], amrex::make_alias, scomp, ncomp))};
        Array<MultiFab,AMREX_SPACEDIM> crse
            {AMREX_D_DECL(MultiFab(*cmf[0], amrex::make_alias, scomp, ncomp),
                          MultiFab(*cmf[1], amrex::make_alias, scomp, ncomp),
                          MultiFab(*cmf[2], amrex::make_alias, scomp, ncomp))};
        amrex::average_down_faces(GetArrOfConstPtrs(fine), GetArrOfPtrs(crse),
                                  IntVect(rr), *cgeom);
    }

    void amrex_fi_average_cellcenter_to_face (MultiFab* fc[], const MultiFab* cc, const Geometry* geom)
    {
        Vector<MultiFab*> fcv {AMREX_D_DECL(fc[0], fc[1], fc[2])};
        amrex::average_cellcenter_to_face(fcv, *cc, *geom);
    }
}
