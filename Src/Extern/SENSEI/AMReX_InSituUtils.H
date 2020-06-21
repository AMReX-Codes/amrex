#ifndef AMReX_InSituUtils_H
#define AMReX_InSituUtils_H

#include <vtkDataSetAttributes.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_MFIter.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>

#include <map>
#include <utility>
#include <string>
#include <vector>

namespace amrex {
namespace InSituUtils {

// The following bit fields are consistent with VisIt ghost zones specification
// For details, see http://www.visitusers.org/index.php?title=Representing_ghost_data

enum CellGhostTypes
{
  DUPLICATECELL           = 1,  // the cell is present on multiple processors
  HIGHCONNECTIVITYCELL    = 2,  // the cell has more neighbors than in a regular mesh
  LOWCONNECTIVITYCELL     = 4,  // the cell has less neighbors than in a regular mesh
  REFINEDCELL             = 8,  // other cells are present that refines it.
  EXTERIORCELL            = 16, // the cell is on the exterior of the data set
  HIDDENCELL              = 32  // the cell is needed to maintain connectivity,
                                // but the data values should be ignored.
};

enum PointGhostTypes
{
  DUPLICATEPOINT          = 1,  // the cell is present on multiple processors
  HIDDENPOINT             = 2   // the point is needed to maintain connectivity,
                                // but the data values should be ignored.
};



// traits helper for mapping between amrex_real and vtkDataArray
template <typename cpp_t> struct amrex_tt {};

#define amrex_tt_specialize(cpp_t, vtk_t, vtk_t_e)      \
template <>                                             \
struct amrex_tt<cpp_t>                                  \
{                                                       \
    using vtk_type = vtk_t;                             \
                                                        \
    static                                              \
    constexpr int vtk_type_enum() { return vtk_t_e; }   \
};

amrex_tt_specialize(float, vtkFloatArray, VTK_FLOAT)
amrex_tt_specialize(double, vtkDoubleArray, VTK_DOUBLE)


// helpers to modify values
template<typename n_t>
struct assignOp { void operator()(n_t &elem, n_t val) { elem = val; } };

template<typename n_t>
struct orOp { void operator()(n_t &elem, n_t val) { elem |= val; } };

// --------------------------------------------------------------------------
template <typename n_t, typename op_t>
int UpdatePatch(long *len, int *lo, int *hi, n_t *darray, n_t val)
{
    op_t updateOp;

    long nx = len[0];
    long nxy = len[0]*len[1];

    for (int k = lo[2]; k <= hi[2]; ++k)
    {
        for (int j = lo[1]; j <= hi[1]; ++j)
        {
            for (int i = lo[0]; i <= hi[0]; ++i)
            {
                updateOp(darray[k*nxy + j*nx + i], val);
            }
        }
    }
    return 0;
}

// --------------------------------------------------------------------------
template <typename n_t>
int AllocateBoxArray(const amrex::Box &pdom, const amrex::BoxArray &boxes,
    const amrex::DistributionMapping &dmap, int ng, std::vector<n_t*> &darrays)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned int nBoxes = boxes.size();

    darrays.resize(nBoxes, nullptr);

    // loop over boxes
    for (unsigned int p = 0; p < nBoxes; ++p)
    {
        // skip non local box
        if (dmap[p] != rank)
            continue;

        // for computing indices
        amrex::Box box = boxes[p];

        // add ghosts
        for (unsigned int q = 0; q < AMREX_SPACEDIM; ++q)
            box.grow(q, ng);

        // TODO -- handle node cenetred data via flag
        // and surroundingNodes

        // array size
        int lo[3] = {AMREX_ARLIM(box.loVect())};
        int hi[3] = {AMREX_ARLIM(box.hiVect())};
        long len[3] = {0};

        for (int q = 0; q < 3; ++q)
            len[q] = hi[q] - lo[q] + 1;

        long nxyz = 1;
        for (int q = 0; q < 3; ++q)
            nxyz *= len[q];

        // allocate array
        n_t *tmp = (n_t*)malloc(sizeof(n_t)*nxyz);

        // default initialize
        for (long q = 0; q < nxyz; ++q)
            tmp[q] = n_t();

        // store at box index
        darrays[p] = tmp;
    }

    return 0;
}

// --------------------------------------------------------------------------
template <typename n_t>
int MaskGhostCells(const amrex::Box &pdom, const amrex::BoxArray &boxes,
    const amrex::DistributionMapping &dmap, unsigned int ng,
    std::vector<n_t*> &darrays)
{
    // no ghosts to mask
    if (ng < 1)
        return 0;

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // loop over boxes
    unsigned int nBoxes = boxes.size();
    for (unsigned int p = 0; p < nBoxes; ++p)
    {
        // skip non local box
        if (dmap[p] != rank)
            continue;

        // for computing indices
        amrex::Box box = boxes[p];

        // add ghosts
        for (int q = 0; q < AMREX_SPACEDIM; ++q)
            box.grow(q, ng);

        // restrict to problem domain
        //box &= pdom;

        // shift to array indexing space
        box.shift(-1*box.smallEnd());

        // array size
        int lo[3] = {AMREX_ARLIM(box.loVect())};
        int hi[3] = {AMREX_ARLIM(box.hiVect())};

        long len[3] = {0};
        for (int q = 0; q < 3; ++q)
            len[q] = hi[q] - lo[q] + 1;

        // get array
        n_t *tmp = darrays[p];

        // for each valid direction
        for (int q = 0; q < AMREX_SPACEDIM; ++q)
        {
            int glo[3] = {0};
            int ghi[3] = {0};

            // fill low side ghosts
            for (int r = 0; r < 3; ++r)
                glo[r] = lo[r];

            for (int r = 0; r < 3; ++r)
                ghi[r] = hi[r];

            ghi[q] = glo[q] + ng - 1;

            UpdatePatch<unsigned char,orOp<unsigned char>>(len,
                glo, ghi, tmp, CellGhostTypes::DUPLICATECELL);

            // fill high side ghosts
            for (int r = 0; r < 3; ++r)
                glo[r] = lo[r];

            for (int r = 0; r < 3; ++r)
                ghi[r] = hi[r];

            glo[q] = ghi[q] + ng - 1;

            UpdatePatch<unsigned char,orOp<unsigned char>>(len,
                glo, ghi, tmp, CellGhostTypes::DUPLICATECELL);
        }
    }

    return 0;
}

// --------------------------------------------------------------------------
template <typename n_t>
int MaskCoveredCells(const amrex::Box &pdom, const amrex::BoxArray &cBoxes,
    const amrex::DistributionMapping &cMap, amrex::BoxArray fBoxes,
    const IntVect &fRefRatio, int ng, std::vector<n_t*> &cMasks)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // work in the coarse level index space
    fBoxes.coarsen(fRefRatio);

    // loop over the local coarse boxes
    unsigned int nCBoxes = cBoxes.size();
    for (unsigned int p = 0; p < nCBoxes; ++p)
    {
        // skip non local coarse box. we do not need to mask it.
        if (cMap[p] != rank)
            continue;

        // this one tells us info for computing indices
        amrex::Box cBox = cBoxes[p];

        // add ghosts
        for (int q = 0; q < AMREX_SPACEDIM; ++q)
            cBox.grow(q, ng);

        // restrict to problem domain
        //cBox &= pdom;

        // array size
        int cLo[3] = {AMREX_ARLIM(cBox.loVect())};
        int cHi[3] = {AMREX_ARLIM(cBox.hiVect())};

        long cLen[3] = {0};
        for (int q = 0; q < 3; ++q)
            cLen[q] = cHi[q] - cLo[q] + 1;

        // mask array
        n_t *cMask = cMasks[p];

        // loop over fine boxes
        unsigned int nFBoxes = fBoxes.size();
        for (unsigned int pp = 0; pp < nFBoxes; ++pp)
        {
            // intersect, the intersection is the covered region
            amrex::Box fBox = fBoxes[pp];;
            amrex::Box iBox = fBox & cBox;


            if (!iBox.isEmpty())
            {
                // shift to array indexing space
                iBox -= cBox.smallEnd();

                // indeces of covered region
                int iLo[3] = {AMREX_ARLIM(iBox.loVect())};
                int iHi[3] = {AMREX_ARLIM(iBox.hiVect())};

                UpdatePatch<n_t, orOp<n_t>>(cLen, iLo, iHi, cMask,
                    CellGhostTypes::REFINEDCELL);
            }
        }
    }

    return 0;
}


// helper to track names and centerings of the avaliable arrays
class StateMap
{
public:
    void Clear(){ this->Map.clear(); }

    int GetIndex(const std::string &arrayName,
       int centering, int &desc, int &comp);

    int Size(int centering){ return this->Map[centering].size(); }

    int GetName(int centering, int id, std::string &name);

protected:
    std::map<int,std::map<std::string,std::pair<int,int>>> Map;
};

}
}

#endif

