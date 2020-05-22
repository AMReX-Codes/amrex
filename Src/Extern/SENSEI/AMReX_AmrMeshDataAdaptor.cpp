#include "AMReX_AmrMeshDataAdaptor.H"

#include "Profiler.h"
#include "Error.h"
#include "VTKUtils.h"

#include <vtkObjectFactory.h>
#include <vtkOverlappingAMR.h>
#include <vtkAMRBox.h>
#include <vtkUniformGrid.h>
#include <vtkDataSetAttributes.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_MFIter.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>
#include <AMReX_InSituUtils.H>

#include <iostream>
#include <sstream>
#include <array>

namespace amrex {
namespace InSituUtils {

// helper to track names and centerings of the avaliable arrays
class MeshStateMap : public amrex::InSituUtils::StateMap
{
public:
    int Initialize(
        const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
        const std::vector<std::vector<std::string>> &names);
};

// --------------------------------------------------------------------------
int MeshStateMap::Initialize(
    const std::vector<amrex::Vector<amrex::MultiFab> *> &states,
    const std::vector<std::vector<std::string>> &names)
{
    int nStates = states.size();
    for (int i = 0; i < nStates; ++i)
    {
        amrex::MultiFab& state = states[i]->at(0);
        int nComp = state.nComp();

        for (int j = 0; j < nComp; ++j)
        {
            const std::string &arrayName = names[i][j];

            if (state.is_cell_centered())
            {
                this->Map[vtkDataObject::CELL][arrayName] = std::make_pair(i,j);
            }
            else if (state.is_nodal())
            {
                this->Map[vtkDataObject::POINT][arrayName] = std::make_pair(i,j);
            }
        }
    }

    return 0;
}
}

// data adaptor's internal data
struct AmrMeshDataAdaptor::InternalsType
{
    InternalsType() : Mesh(nullptr), PinMesh(0),
        States(), StateMetadata() {}

    amrex::AmrMesh *Mesh;
    int PinMesh;
    std::vector<amrex::Vector<amrex::MultiFab> *> States;
    std::vector<std::vector<std::string>> Names;
    amrex::InSituUtils::MeshStateMap StateMetadata;
#if SENSEI_VERSION_MAJOR < 3
    std::vector<vtkDataObject*> ManagedObjects;
#endif
};

//-----------------------------------------------------------------------------
senseiNewMacro(AmrMeshDataAdaptor);

//-----------------------------------------------------------------------------
AmrMeshDataAdaptor::AmrMeshDataAdaptor() :
    Internals(new AmrMeshDataAdaptor::InternalsType())
{
}

//-----------------------------------------------------------------------------
AmrMeshDataAdaptor::~AmrMeshDataAdaptor()
{
    delete this->Internals;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::SetDataSource(amrex::AmrMesh *mesh,
    const std::vector<amrex::Vector<amrex::MultiFab>*> &states,
    const std::vector<std::vector<std::string>> &names)
{
    this->ReleaseData();

    this->Internals->Mesh = mesh;
    this->Internals->States = states;
    this->Internals->Names = names;

    this->Internals->StateMetadata.Initialize(states, names);

    return 0;
}

//-----------------------------------------------------------------------------
void AmrMeshDataAdaptor::SetPinMesh(int pinMesh)
{
    this->Internals->PinMesh = pinMesh;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes)
{
    numMeshes = 1;
    return 0;
}



#if SENSEI_VERSION_MAJOR >= 3
int AmrMeshDataAdaptor::GetMeshMetadata(unsigned int id,
    sensei::MeshMetadataPtr &metadata)
{
    sensei::TimeEvent<64> event("AmrMeshDataAdaptor::GetMeshMetadata");

    if (id != 0)
      {
      SENSEI_ERROR("invalid mesh id " << id)
      return -1;
      }

    // AMR data is always expected to be a global view
    metadata->GlobalView = true;

    metadata->MeshName = "mesh";
    metadata->MeshType = VTK_OVERLAPPING_AMR;
    metadata->BlockType = VTK_UNIFORM_GRID;
    metadata->NumBlocks = 0;
    metadata->NumCells = 0;
    metadata->NumPoints = 0;
    metadata->NumBlocksLocal = {-1};
    metadata->CoordinateType = InSituUtils::amrex_tt<amrex_real>::vtk_type_enum();
    metadata->StaticMesh = 0;

    // num levels
    metadata->NumLevels = this->Internals->Mesh->finestLevel() + 1;

    // total number of blocks
    for (unsigned int i = 0; i < metadata->NumLevels; ++i)
    {
        int numBlocks = this->Internals->Mesh->boxArray(i).size();
        metadata->NumBlocks += numBlocks;
    }

    // allocate per block decomp metadata
    if (metadata->Flags.BlockDecompSet())
    {
        metadata->BlockIds.resize(metadata->NumBlocks);
        metadata->BlockOwner.resize(metadata->NumBlocks);
        metadata->BlockLevel.resize(metadata->NumBlocks);
        metadata->RefRatio.resize(metadata->NumLevels);
        metadata->BlocksPerLevel.resize(metadata->NumLevels);
    }

    // allocate per block size metadata
    if (metadata->Flags.BlockSizeSet())
    {
        metadata->NumCells = 0;
        metadata->NumPoints = 0;

        metadata->BlockNumCells.resize(metadata->NumBlocks);
        metadata->BlockNumPoints.resize(metadata->NumBlocks);
    }

    // allocate per block bounds metadata
    if (metadata->Flags.BlockBoundsSet())
    {
        metadata->BlockBounds.resize(metadata->NumBlocks);
    }

    // allocate per block extent metadata
    if (metadata->Flags.BlockBoundsSet())
    {
        metadata->BlockExtents.resize(metadata->NumBlocks);
    }

    // global bounds
    const amrex::RealBox& pd = this->Internals->Mesh->Geom(0).ProbDomain();

    double pdLo[3] = {AMREX_ARLIM(pd.lo())};
    double pdHi[3] = {AMREX_ARLIM(pd.hi())};

    // PinMesh works around a bug in VisIt 2.13.2.
    // force the origin to 0,0,0
    if (this->Internals->PinMesh)
    {
        for (int i = 0; i < 3; ++i)
        {
            pdHi[i] = pdHi[i] - pdLo[i];
            pdLo[i] = 0.0;
        }
    }

    // global bounds
    if (metadata->Flags.BlockBoundsSet())
    {
        metadata->Bounds = std::array<double,6>(
            {pdLo[0], pdHi[0], pdLo[1], pdHi[1], pdLo[2], pdHi[2]});
    }

    // global extent (note: VTK uses point centered indexing)
    const amrex::Box& cdom = this->Internals->Mesh->Geom(0).Domain();
    amrex::Box ndom = surroundingNodes(cdom);

    int i0[3] = {AMREX_ARLIM(ndom.loVect())};
    int i1[3] = {AMREX_ARLIM(ndom.hiVect())};

    // global extent
    if (metadata->Flags.BlockExtentsSet())
    {
        metadata->Extent = std::array<int,6>(
            {i0[0], i1[0], i0[1], i1[1], i0[2], i1[2]});
    }

    // num ghosts
    metadata->NumGhostCells = this->Internals->States[0]->at(0).nGrow();
    metadata->NumGhostNodes = 0;

    // per array metadata
    amrex::MultiFab& state0 = this->Internals->States[id]->at(0);

    // number of arrays
    metadata->NumArrays = state0.nComp();

    // resize per array containers
    metadata->ArrayName.resize(metadata->NumArrays);
    metadata->ArrayComponents.resize(metadata->NumArrays);
    metadata->ArrayType.resize(metadata->NumArrays);
    metadata->ArrayCentering.resize(metadata->NumArrays);

    for (int j = 0; j < metadata->NumArrays; ++j)
    {
        // array name
        metadata->ArrayName[j] = this->Internals->Names[id][j];
        // scalar, vector, tensor
        metadata->ArrayComponents[j] = 1;
        // POD type
        metadata->ArrayType[j] = InSituUtils::amrex_tt<amrex_real>::vtk_type_enum();
        // mesh centering
        if (state0.is_cell_centered())
        {
            metadata->ArrayCentering[j] = vtkDataObject::CELL;
        }
        else if (state0.is_nodal())
        {
            metadata->ArrayCentering[j] = vtkDataObject::POINT;
        }
        else
        {
            metadata->ArrayCentering[j] = vtkDataObject::FIELD;
        }
    }

    // TODO -- port code from AmrDataAdaptor
    if (metadata->Flags.BlockArrayRangeSet())
    {
        SENSEI_ERROR("Per-array per-block array ranges not yet implemented")
    }

    // gather per block metadata
    long gid = 0;
    for (unsigned int i = 0; i < metadata->NumLevels; ++i)
    {
        // domain decomp
        const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

        // ghost zones
        amrex::MultiFab& state = this->Internals->States[0]->at(i);
        unsigned int ng = state.nGrow();

        // spacing
        const amrex::Geometry &geom = this->Internals->Mesh->Geom(i);
        double spacing [3] = {AMREX_ARLIM(geom.CellSize())};

        // boxes
        const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
        unsigned int nBoxes = ba.size();

        if (metadata->Flags.BlockDecompSet())
        {
            // blocks per level
            metadata->BlocksPerLevel[i] = nBoxes;

            // refinement ratio
            if (i != this->Internals->Mesh->maxLevel())
                metadata->RefRatio[i] = std::array<int,3>(
                    {AMREX_ARLIM(this->Internals->Mesh->refRatio(i))});
            else
                metadata->RefRatio[i] = std::array<int,3>({1,1,1});
        }

        // process boxes this level
        for (unsigned int j = 0; j < nBoxes; ++j, ++gid)
        {
            // cell centered box
            amrex::Box cbox = ba[j];

            // add ghost zones
            for (int q = 0; q < AMREX_SPACEDIM; ++q)
                cbox.grow(q, ng);

            // node centered box
            amrex::Box nbox = surroundingNodes(cbox);

            // node centered dimensions
            int nboxLo[3] = {AMREX_ARLIM(nbox.loVect())};
            int nboxHi[3] = {AMREX_ARLIM(nbox.hiVect())};

            // domain decomp
            if (metadata->Flags.BlockDecompSet())
            {
                metadata->BlockOwner[gid] = dmap[j];
                metadata->BlockIds[gid] = gid;

                metadata->BlockLevel[gid] = i;
            }

            // block sizes
            if (metadata->Flags.BlockSizeSet())
            {
                long nPts = nbox.numPts();
                long nCells = cbox.numPts();

                metadata->NumPoints += nPts;
                metadata->NumCells += nCells;

                metadata->BlockNumPoints[gid] = nPts;
                metadata->BlockNumCells[gid] = nCells;
            }

            // block extent
            if (metadata->Flags.BlockExtentsSet())
            {
                metadata->BlockExtents[gid] = std::array<int,6>({
                    nboxLo[0], nboxHi[0], nboxLo[1], nboxHi[1], nboxLo[2], nboxHi[2]});
            }

            // block bounds
            if (metadata->Flags.BlockBoundsSet())
            {
                metadata->BlockBounds[gid] = std::array<double,6>(
                    {pdLo[0] + spacing[0]*nboxLo[0], pdHi[0] + spacing[0]*nboxHi[0],
                    pdLo[1] + spacing[1]*nboxLo[1], pdHi[1] + spacing[1]*nboxHi[1],
                    pdLo[2] + spacing[2]*nboxLo[2], pdHi[2] + spacing[2]*nboxHi[2]});
            }

        }
    }

    return 0;
}
#else

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetMeshName(unsigned int id, std::string &meshName)
{
    meshName = "mesh";
    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetNumberOfArrays(const std::string &meshName,
    int association, unsigned int &numberOfArrays)
{
    numberOfArrays = 0;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    if ((association != vtkDataObject::POINT) &&
        (association != vtkDataObject::CELL))
    {
        SENSEI_ERROR("Invalid association " << association)
        return -1;
    }

    if (this->Internals->States.size() < 1)
    {
        SENSEI_ERROR("No simualtion data, missing call to SetDataSource?")
        return -1;
    }

    numberOfArrays = this->Internals->StateMetadata.Size(association);

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetArrayName(const std::string &meshName,
    int association, unsigned int index, std::string &arrayName)
{
    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    if (this->Internals->StateMetadata.GetName(association, index, arrayName))
    {
        SENSEI_ERROR("No array named \"" << arrayName << "\" in "
            << sensei::VTKUtils::GetAttributesName(association)
            << " data")
        return -1;
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetMeshHasGhostNodes(const std::string &meshName, int &nLayers)
{
    nLayers = 0;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetMeshHasGhostCells(const std::string &meshName, int &nLayers)
{
    nLayers = 0;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    if (this->Internals->States.size() < 1)
    {
        SENSEI_ERROR("No simualtion data, missing call to SetDataSource?")
        return -1;
    }

    nLayers = 1;

    return 0;
}
#endif


//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::GetMesh(const std::string &meshName,
    bool structureOnly, vtkDataObject *&mesh)
{
    mesh = nullptr;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("No mesh named \"" << meshName << "\"")
        return -1;
    }

    if (this->Internals->States.size() < 1)
    {
        SENSEI_ERROR("No simualtion data, missing call to SetDataSource?")
        return -1;
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

    // initialize new vtk datasets
    vtkOverlappingAMR *amrMesh = vtkOverlappingAMR::New();
#if SENSEI_VERSION_MAJOR < 3
    Internals->ManagedObjects.push_back(amrMesh);
#endif
    mesh = amrMesh;

    // num levels and blocks per level
    std::vector<int> nBlocks(nLevels);
    for (unsigned int i = 0; i < nLevels; ++i)
        nBlocks[i] = this->Internals->Mesh->boxArray(i).size();

    amrMesh->Initialize(nLevels, nBlocks.data());

    // origin
    const amrex::RealBox& pd = this->Internals->Mesh->Geom(0).ProbDomain();
    double origin[3] = {AMREX_ARLIM(pd.lo())};

    // PinMesh works around a bug in VisIt 2.13.2.
    // force the origin to 0,0,0
    if (this->Internals->PinMesh)
    {
        for (int i = 0; i < 3; ++i)
            origin[i] = 0.0;
    }

    amrMesh->SetOrigin(origin);

    long gid = 0;
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        // domain decomp
        const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

        // ghost zones
        amrex::MultiFab& state = this->Internals->States[0]->at(i);
        unsigned int ng = state.nGrow();

        // spacing
        const amrex::Geometry &geom = this->Internals->Mesh->Geom(i);
        double spacing [3] = {AMREX_ARLIM(geom.CellSize())};
        amrMesh->SetSpacing(i, spacing);

        // refinement ratio
        int rr = 1;
        if (i != this->Internals->Mesh->maxLevel())
            rr = this->Internals->Mesh->refRatio(i)[0];
        amrMesh->SetRefinementRatio(i, rr);

        // loop over boxes
        const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
        unsigned int nBoxes = ba.size();

        for (unsigned int j = 0; j < nBoxes; ++j)
        {
            // cell centered box
            amrex::Box cbox = ba[j];

            // cell centered dimensions
            int cboxLo[3] = {AMREX_ARLIM(cbox.loVect())};
            int cboxHi[3] = {AMREX_ARLIM(cbox.hiVect())};

            // vtk's representation of box metadata
            vtkAMRBox block(cboxLo, cboxHi);
            amrMesh->SetAMRBox(i, j, block);
            amrMesh->SetAMRBlockSourceIndex(i, j, gid++);

            // skip building a vtk amrMesh for the non local boxes
            if (dmap[j] != rank)
                continue;

            // add ghost zones
            for (int q = 0; q < AMREX_SPACEDIM; ++q)
                cbox.grow(q, ng);

            // node centered box
            amrex::Box nbox = surroundingNodes(cbox);

            // node centered dimensions
            int nboxLo[3] = {AMREX_ARLIM(nbox.loVect())};
            int nboxHi[3] = {AMREX_ARLIM(nbox.hiVect())};

            // new vtk uniform amrMesh, node centered
            vtkUniformGrid *ug = vtkUniformGrid::New();
            ug->SetOrigin(origin);
            ug->SetSpacing(spacing);
            ug->SetExtent(nboxLo[0], nboxHi[0],
                nboxLo[1], nboxHi[1], nboxLo[2], nboxHi[2]);

            // pass the block into vtk
            amrMesh->SetDataSet(i, j, ug);
            ug->Delete();
        }
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::AddGhostNodesArray(vtkDataObject *mesh,
    const std::string &meshName)
{
    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::AddGhostCellsArray(vtkDataObject* mesh,
    const std::string &meshName)
{
    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    vtkOverlappingAMR *amrMesh = dynamic_cast<vtkOverlappingAMR*>(mesh);
    if (!amrMesh)
    {
        SENSEI_ERROR("Invalid mesh type "
            << (mesh ? mesh->GetClassName() : "nullptr"))
    }

    if (this->Internals->States.size() < 1)
    {
        SENSEI_ERROR("No simualtion data, missing call to SetDataSource?")
        return -1;
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // loop over levels
    unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

    std::vector<std::vector<unsigned char*>> masks(nLevels);
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        // allocate mask arrays
        const amrex::BoxArray &boxes = this->Internals->Mesh->boxArray(i);
        const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);
        const amrex::Box &pdom = this->Internals->Mesh->Geom(i).Domain();

        amrex::MultiFab& state = this->Internals->States[0]->at(i);
        unsigned int ng = state.nGrow();

        std::vector<unsigned char*> mask;
        InSituUtils::AllocateBoxArray<unsigned char>(
            pdom, boxes, dmap, ng, mask);

        // mask ghost cells
        InSituUtils::MaskGhostCells<unsigned char>(
            pdom, boxes, dmap, ng, mask);

        // store mask array
        masks[i] = mask;
    }

    // loop over coarse levels
    unsigned int nCoarseLevels = nLevels - 1;
    for (unsigned int i = 0; i < nCoarseLevels; ++i)
    {
        int ii = i + 1;

        // mask regions covered by refinement
        amrex::MultiFab& state = this->Internals->States[0]->at(i);
        unsigned int ng = state.nGrow();

        const amrex::Box &pdom = this->Internals->Mesh->Geom(i).Domain();
        const amrex::BoxArray &cBoxes = this->Internals->Mesh->boxArray(i);
        const amrex::DistributionMapping &cMap = this->Internals->Mesh->DistributionMap(i);
        const amrex::BoxArray &fBoxes = this->Internals->Mesh->boxArray(ii);
        amrex::IntVect cRefRatio = this->Internals->Mesh->refRatio(i);

        InSituUtils::MaskCoveredCells<unsigned char>(
            pdom, cBoxes, cMap, fBoxes, cRefRatio, ng, masks[i]);
    }

    // loop over levels
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

        // mask arrays for this level
        std::vector<unsigned char*> &mask = masks[i];

        // loop over boxes
        const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
        unsigned int nBoxes = ba.size();

        for (unsigned int j = 0; j < nBoxes; ++j)
        {
            if (dmap[j] != rank)
                continue;

            vtkUniformGrid *blockMesh = amrMesh->GetDataSet(i, j);

            if (!blockMesh)
            {
                SENSEI_ERROR("Empty block " << i << ", " << j)
                return -1;
            }

            long nCells = blockMesh->GetNumberOfCells();

            // transfer mask array into vtk
            vtkUnsignedCharArray *ga = vtkUnsignedCharArray::New();
            ga->SetName("vtkGhostType");
            ga->SetArray(mask[j], nCells, 0);
            blockMesh->GetCellData()->AddArray(ga);
            ga->Delete();

            // for debug can visualize the ghost cells
            // FIXME -- a bug in Catalyst ignores internal ghost zones
            // when using the VTK writrer. Until that bug gets fixed, one
            // can manually inject this copy using a PV Python filter
            ga = vtkUnsignedCharArray::New();
            ga->SetName("GhostType");
            ga->SetArray(mask[j], nCells, 1);
            blockMesh->GetCellData()->AddArray(ga);
            ga->Delete();
        }
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::AddArray(vtkDataObject* mesh,
    const std::string &meshName, int association,
    const std::string &arrayName)
{
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    vtkOverlappingAMR *amrMesh = dynamic_cast<vtkOverlappingAMR*>(mesh);
    if (!amrMesh)
    {
        SENSEI_ERROR("Invalid mesh type "
            << (mesh ? mesh->GetClassName() : "nullptr"))
    }

    if (this->Internals->States.size() < 1)
    {
        SENSEI_ERROR("No simualtion data, missing call to SetDataSource?")
        return -1;
    }

    if ((association != vtkDataObject::CELL) &&
        (association != vtkDataObject::CELL))
    {
        SENSEI_ERROR("Invalid association " << association)
        return -1;
    }

    // find the indices of the multifab and component within for
    // the named array
    int fab = 0;
    int comp = 0;
    if (this->Internals->StateMetadata.GetIndex(arrayName, association, fab, comp))
    {
        SENSEI_ERROR("Failed to locate descriptor for "
            << sensei::VTKUtils::GetAttributesName(association)
            << " data array \"" << arrayName << "\"")
        return -1;
    }

    // loop over levels
    unsigned int nLevels = this->Internals->Mesh->finestLevel() + 1;

    for (unsigned int i = 0; i < nLevels; ++i)
    {
        // domain decomp
        const amrex::DistributionMapping &dmap = this->Internals->Mesh->DistributionMap(i);

        // ghost zones
        amrex::MultiFab& state = this->Internals->States[fab]->at(i);
        unsigned int ng = state.nGrow();

        // check centering
        if (!((association == vtkDataObject::CELL) && state.is_cell_centered()) &&
            !((association == vtkDataObject::POINT) && state.is_nodal()))
        {
            SENSEI_ERROR("association does not match MultiFab centering")
            return -1;
        }

        // check component id
        int nComp = state.nComp();
        if (comp >= nComp)
        {
            SENSEI_ERROR("Component " << comp << " out of bounds")
            return -1;
        }

        // loop over boxes
        const amrex::BoxArray& ba = this->Internals->Mesh->boxArray(i);
        unsigned int nBoxes = ba.size();

        for (unsigned int j = 0; j < nBoxes; ++j)
        {
            // cell centered box
            amrex::Box cbox = ba[j];

            // add ghost zones
            for (int q = 0; q < AMREX_SPACEDIM; ++q)
                cbox.grow(q, ng);

            // cell centered dimensions
            int cboxLo[3] = {AMREX_ARLIM(cbox.loVect())};
            int cboxHi[3] = {AMREX_ARLIM(cbox.hiVect())};

            // skip building a vtk mesh for the non local boxes
            if (dmap[j] != rank)
                continue;

            // node centered box
            amrex::Box nbox = surroundingNodes(cbox);

            // node centered dimensions
            int nboxLo[3] = {AMREX_ARLIM(nbox.loVect())};
            int nboxHi[3] = {AMREX_ARLIM(nbox.hiVect())};

            // get the block mesh
            vtkUniformGrid *ug = amrMesh->GetDataSet(i, j);

            // node centered size
            long nlen = 1;
            for (int p = 0; p < 3; ++p)
                nlen *= nboxHi[p] - nboxLo[p] + 1;

            // cell centered size
            long clen = 1;
            for (int p = 0; p < 3; ++p)
                clen *= cboxHi[p] - cboxLo[p] + 1;

            // pointer to the data
            amrex_real *pcd = state[j].dataPtr(comp);

            // allocate vtk array
            InSituUtils::amrex_tt<amrex_real>::vtk_type *da =
                InSituUtils::amrex_tt<amrex_real>::vtk_type::New();

            // set component name
            da->SetName(arrayName.c_str());

            if (state[j].box().ixType() == amrex::IndexType::TheCellType())
            {
                // zero copy cell centered
                da->SetArray(pcd, clen, 1);
                ug->GetCellData()->AddArray(da);
            }
            else if (state[j].box().ixType() == amrex::IndexType::TheNodeType())
            {
                // zero copy point centered
                da->SetArray(pcd, nlen, 1);
                ug->GetPointData()->AddArray(da);
            }
            else
            {
                SENSEI_WARNING("Face or edge centered component " << comp << " skipped")
            }

            da->Delete();

#if defined(SENSEI_DEBUG)
            // mark level id
            vtkFloatArray *la = vtkFloatArray::New();
            la->SetName("amrex_level_id");
            la->SetNumberOfTuples(clen);
            la->Fill(i);
            ug->GetCellData()->AddArray(la);
            la->Delete();

            // mark mpi rank
            vtkFloatArray *ra = vtkFloatArray::New();
            ra->SetName("amrex_mpi_rank");
            ra->SetNumberOfTuples(clen);
            ra->Fill(rank);
            ug->GetCellData()->AddArray(ra);
            ra->Delete();
#endif
        }
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrMeshDataAdaptor::ReleaseData()
{
    this->Internals->Mesh = nullptr;
    this->Internals->States.clear();
    this->Internals->Names.clear();
    this->Internals->StateMetadata.Clear();

#if SENSEI_VERSION_MAJOR < 3
    // free up mesh objects we allocated
    size_t n = this->Internals->ManagedObjects.size();
    for (size_t i = 0; i < n; ++i)
        this->Internals->ManagedObjects[i]->Delete();
    this->Internals->ManagedObjects.clear();
#endif

    return 0;
}

}
