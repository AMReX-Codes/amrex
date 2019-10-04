#include "AMReX_AmrDataAdaptor.H"

#include "MPIUtils.h"
#include "STLUtils.h"
#include "VTKUtils.h"
#include "Profiler.h"
#include "Error.h"

#include <vtkObjectFactory.h>
#include <vtkOverlappingAMR.h>
#include <vtkAMRBox.h>
#include <vtkUniformGrid.h>
#include <vtkXMLUniformGridAMRWriter.h>
#include <vtkDataSetAttributes.h>
#include <vtkUnsignedCharArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include <vtkCellData.h>
#include <vtkPointData.h>

#include <AMReX_AmrLevel.H>
#include <AMReX_BoxArray.H>
#include <AMReX_Geometry.H>
#include <AMReX_RealBox.H>
#include <AMReX_Box.H>
#include <AMReX_ArrayLim.H>
#include <AMReX_StateData.H>
#include <AMReX_MFIter.H>
#include <AMReX_DistributionMapping.H>
#include <AMReX_IndexType.H>
#include <AMReX_InSituUtils.H>

#include <iostream>
#include <sstream>
#include <map>
#include <utility>

// return the number of levels currently in use
static
unsigned int numActiveLevels(
    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels)
{
    unsigned int nLevels = levels.size();
    for (int i = 0; i < nLevels; ++i)
    {
        if (!levels[i])
        {
            nLevels = i;
            break;
        }
    }
    return nLevels;
}


namespace amrex {
namespace InSituUtils {

// helper to track names and centerings of the avaliable arrays
class DescriptorMap : public amrex::InSituUtils::StateMap
{
public:
    int Initialize(const DescriptorList &descriptors);
};

// --------------------------------------------------------------------------
int DescriptorMap::Initialize(const DescriptorList &descriptors)
{
    int ndesc = descriptors.size();
    for (int i = 0; i < ndesc; ++i)
    {
        const StateDescriptor &desc = descriptors[i];

        int ncomp = desc.nComp();
        IndexType itype = desc.getType();

        for (int j = 0; j < ncomp; ++j)
        {
            std::string arrayName = desc.name(j);

            if (itype.cellCentered())
            {
                this->Map[vtkDataObject::CELL][arrayName] = std::make_pair(i,j);
            }
            else if (itype.nodeCentered())
            {
                this->Map[vtkDataObject::POINT][arrayName] = std::make_pair(i,j);
            }
        }
    }

    return 0;
}
}

// data adaptor's internal data
struct AmrDataAdaptor::InternalsType
{
    InternalsType() : SimData(nullptr), PinMesh(0) {}

    amrex::Amr *SimData;
    int PinMesh;
    amrex::InSituUtils::DescriptorMap SimMetadata;
#if SENSEI_VERSION_MAJOR < 3
    std::vector<vtkDataObject*> ManagedObjects;
#endif
};

//-----------------------------------------------------------------------------
senseiNewMacro(AmrDataAdaptor);

//-----------------------------------------------------------------------------
AmrDataAdaptor::AmrDataAdaptor() :
    Internals(new AmrDataAdaptor::InternalsType())
{
}

//-----------------------------------------------------------------------------
AmrDataAdaptor::~AmrDataAdaptor()
{
    delete this->Internals;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::SetDataSource(amrex::Amr *amr)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::SetDataSource");

    this->ReleaseData();

    this->Internals->SimData = amr;

    // array metadata
    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
        this->Internals->SimData->getAmrLevels();

    this->Internals->SimMetadata.Initialize(levels[0]->get_desc_lst());

    return 0;
}

//-----------------------------------------------------------------------------
void AmrDataAdaptor::SetPinMesh(int pinMesh)
{
    this->Internals->PinMesh = pinMesh;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetNumberOfMeshes(unsigned int &numMeshes)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::GetNumberOfMeshes");
    numMeshes = 1;
    return 0;
}

#if SENSEI_VERSION_MAJOR >= 3
//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetMeshMetadata(unsigned int id,
    sensei::MeshMetadataPtr &metadata)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::GetMeshMetadata");

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
    metadata->NumBlocksLocal = {-1};
    metadata->CoordinateType = InSituUtils::amrex_tt<amrex_real>::vtk_type_enum();
    metadata->StaticMesh = 0;

    // TODO
    //metadata->PeriodicBoundary = ;

    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
     this->Internals->SimData->getAmrLevels();

    // num levels and blocks per level
    metadata->NumLevels = numActiveLevels(levels);

    metadata->NumBlocks = 0;
    metadata->BlocksPerLevel.resize(metadata->NumLevels);
    for (unsigned int i = 0; i < metadata->NumLevels; ++i)
    {
        unsigned long nb = levels[i]->boxArray().size();
        metadata->NumBlocks += nb;
        metadata->BlocksPerLevel[i] = nb;
    }

    // bounds
    const amrex::RealBox& pd = levels[0]->Geom().ProbDomain();

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

    if (metadata->Flags.BlockBoundsSet())
    {
        metadata->Bounds = std::array<double,6>({pdLo[0], pdHi[0],
            pdLo[1], pdHi[1], pdLo[2], pdHi[2]});

        metadata->BlockBounds.reserve(metadata->NumBlocks);
    }

    // extent
    const amrex::Box& dom = levels[0]->Geom().Domain();
    int domLo[3] = {AMREX_ARLIM(dom.smallEnd())};
    int domHi[3] = {AMREX_ARLIM(dom.bigEnd())};

    if (metadata->Flags.BlockExtentsSet())
    {
        metadata->Extent = std::array<int,6>({domLo[0], domHi[0],
            domLo[1], domHi[1], domLo[2], domHi[2]});

        metadata->BlockExtents.reserve(metadata->NumBlocks);
        metadata->BlockLevel.reserve(metadata->NumBlocks);
    }

    // ghost zones
    metadata->NumGhostCells = levels[0]->get_new_data(0).nGrow();

    // arrays
    metadata->NumArrays = 0;
    const DescriptorList &descriptors = levels[0]->get_desc_lst();
    int ndesc = descriptors.size();
    for (int i = 0; i < ndesc; ++i)
    {
        const StateDescriptor &desc = descriptors[i];

        int ncomp = desc.nComp();
        metadata->NumArrays += ncomp;

        IndexType itype = desc.getType();

        for (int j = 0; j < ncomp; ++j)
        {
            std::string arrayName = desc.name(j);
            metadata->ArrayName.push_back(arrayName);
            metadata->ArrayComponents.push_back(1);
            metadata->ArrayType.push_back(InSituUtils::amrex_tt<amrex_real>::vtk_type_enum());

            if (itype.cellCentered())
                metadata->ArrayCentering.push_back(vtkDataObject::CELL);
            else if (itype.nodeCentered())
                metadata->ArrayCentering.push_back(vtkDataObject::POINT);
            else
                metadata->ArrayCentering.push_back(vtkDataObject::FIELD);
        }

    }

    if (metadata->Flags.BlockArrayRangeSet())
        metadata->BlockArrayRange.reserve(metadata->NumBlocks);

    int rank = 0;
    MPI_Comm_rank(this->GetCommunicator(), &rank);

    // per-level and per-block metadata
    long gid = 0;
    for (unsigned int i = 0; i < metadata->NumLevels; ++i)
    {
        // domain decomp
        const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();

        // ghost zones
        unsigned int ng = levels[i]->get_new_data(0).nGrow();

        // spacing
        const amrex::Geometry &geom = levels[i]->Geom();
        double spacing[3] = {AMREX_ARLIM(geom.CellSize())};

        // refinement ratio
        std::array<int,3> rr = {AMREX_ARLIM(levels[i]->fineRatio())};
        metadata->RefRatio.push_back(rr);

        // loop over boxes
        const amrex::BoxArray& ba = levels[i]->boxArray();
        unsigned int nBoxes = ba.size();

        for (unsigned int j = 0; j < nBoxes; ++j)
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
                metadata->BlockOwner.push_back(dmap[j]);
                metadata->BlockIds.push_back(gid++);
            }

            // block sizes
            if (metadata->Flags.BlockSizeSet())
            {
                metadata->BlockNumPoints.push_back(nbox.numPts());
                metadata->BlockNumCells.push_back(cbox.numPts());
            }

            // block extent
            if (metadata->Flags.BlockExtentsSet())
            {
                metadata->BlockExtents.push_back({nboxLo[0], nboxHi[0],
                    nboxLo[1], nboxHi[1], nboxLo[2], nboxHi[2]});

                metadata->BlockLevel.push_back(i);
            }

            // block bounds
            if (metadata->Flags.BlockBoundsSet())
                metadata->BlockBounds.push_back({pdLo[0] + spacing[0]*nboxLo[0],
                    pdLo[0] + spacing[0]*nboxHi[0], pdLo[1] + spacing[1]*nboxLo[1],
                    pdLo[1] + spacing[1]*nboxHi[1], pdLo[2] + spacing[2]*nboxLo[2],
                    pdLo[2] + spacing[2]*nboxHi[2]});

            // only for local blocks
            if ((dmap[j] == rank) && (metadata->Flags.BlockArrayRangeSet()))
            {
                std::vector<std::array<double,2>> arrayRange;
                arrayRange.reserve(metadata->NumArrays);

                // block array range
                for (int k = 0; k < ndesc; ++k)
                {
                    const StateDescriptor &desc = descriptors[k];
                    amrex::MultiFab &state = levels[i]->get_new_data(k);

                    int ncomp = desc.nComp();
                    IndexType itype = desc.getType();

                    for (int l = 0; l < ncomp; ++l)
                    {
                        // calculate min/max on this block for this array
                        amrex_real mn = state[j].min(l);
                        amrex_real mx = state[j].max(l);
                        arrayRange.push_back({mn, mx});
                    }
                }

                metadata->BlockArrayRange.push_back(arrayRange);
            }
        }
    }

    // make the block array range global
    if (metadata->Flags.BlockArrayRangeSet())
    {
        sensei::MPIUtils::GlobalViewV(this->GetCommunicator(), metadata->BlockArrayRange);
        sensei::STLUtils::ReduceRange(metadata->BlockArrayRange, metadata->ArrayRange);
    }

    return 0;
}

#else

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetMeshHasGhostNodes(const std::string &meshName, int &nLayers)
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
int AmrDataAdaptor::GetMeshHasGhostCells(const std::string &meshName, int &nLayers)
{
    nLayers = 0;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    if (!this->Internals->SimData)
    {
        SENSEI_ERROR("No simulation data")
        return -1;
    }

    nLayers = 1;

    return 0;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetNumberOfArrays(const std::string &meshName,
    int association, unsigned int &numberOfArrays)
{
    timer::MarkEvent("AmrDataAdaptor::GetNumberOfArrays");

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

    if (!this->Internals->SimData)
    {
        SENSEI_ERROR("No simulation data")
        return -1;
    }

    numberOfArrays = this->Internals->SimMetadata.Size(association);

    return 0;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetArrayName(const std::string &meshName,
    int association, unsigned int index, std::string &arrayName)
{
    timer::MarkEvent("AmrDataAdaptor::GetArrayName");

    if (meshName != "mesh")
    {
        SENSEI_ERROR("No mesh named \"" << meshName << "\"")
        return -1;
    }

    if (this->Internals->SimMetadata.GetName(association, index, arrayName))
    {
        SENSEI_ERROR("No array named \"" << arrayName << "\" in "
            << sensei::VTKUtils::GetAttributesName(association)
            << " data")
        return -1;
    }

    return 0;
}


#endif

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetMesh(const std::string &meshName,
    bool structureOnly, vtkDataObject *&mesh)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::GetMesh");

    mesh = nullptr;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    int rank = 0;
    MPI_Comm_rank(this->GetCommunicator(), &rank);

    // get levels
    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
     this->Internals->SimData->getAmrLevels();

    unsigned int nLevels = numActiveLevels(levels);

    // initialize new vtk datasets
    vtkOverlappingAMR *amrMesh = vtkOverlappingAMR::New();
#if SENSEI_VERSION_MAJOR < 3
    Internals->ManagedObjects.push_back(amrMesh);
#endif
    mesh = amrMesh;

    // num levels and blocks per level
    std::vector<int> nBlocks(nLevels);
    for (unsigned int i = 0; i < nLevels; ++i)
        nBlocks[i] = levels[i]->boxArray().size();

    amrMesh->Initialize(nLevels, nBlocks.data());

    // origin
    const amrex::RealBox& pd = levels[0]->Geom().ProbDomain();
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
        const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();

        // ghost zones
        amrex::MultiFab &state = levels[i]->get_new_data(0);
        unsigned int ng = state.nGrow();

        // spacing
        const amrex::Geometry &geom = levels[i]->Geom();
        double spacing [3] = {AMREX_ARLIM(geom.CellSize())};
        amrMesh->SetSpacing(i, spacing);

        // refinement ratio
        amrMesh->SetRefinementRatio(i, levels[i]->fineRatio()[0]);

        // loop over boxes
        const amrex::BoxArray& ba = levels[i]->boxArray();
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
int AmrDataAdaptor::AddGhostCellsArray(vtkDataObject* mesh,
    const std::string &meshName)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::AddGhostCellsArray");

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

    if (!this->Internals->SimData)
    {
        SENSEI_ERROR("No simulation data")
        return -1;
    }

    int rank = 0;
    MPI_Comm_rank(this->GetCommunicator(), &rank);

    // loop over levels
    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
     this->Internals->SimData->getAmrLevels();

    unsigned int nLevels = numActiveLevels(levels);

    std::vector<std::vector<unsigned char*>> masks(nLevels);
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        // allocate mask arrays
        const amrex::BoxArray &boxes = levels[i]->boxArray();
        const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();
        const amrex::Box &pdom = levels[i]->Domain();

        amrex::MultiFab& state = levels[i]->get_new_data(0);
        unsigned int ng = state.nGrow();

        std::vector<unsigned char*> mask;
        InSituUtils::AllocateBoxArray<unsigned char>(pdom, boxes, dmap, ng, mask);

        // mask ghost cells
        InSituUtils::MaskGhostCells<unsigned char>(pdom, boxes, dmap, ng, mask);

        // store mask array
        masks[i] = mask;
    }

    // loop over coarse levels
    unsigned int nCoarseLevels = nLevels - 1;
    for (unsigned int i = 0; i < nCoarseLevels; ++i)
    {
        int ii = i + 1;

        // mask regions covered by refinement
        amrex::MultiFab& state = levels[i]->get_new_data(0);
        unsigned int ng = state.nGrow();

        const amrex::Box &pdom = levels[i]->Domain();
        const amrex::BoxArray &cBoxes = levels[i]->boxArray();
        const amrex::DistributionMapping &cMap = levels[i]->DistributionMap();
        const amrex::BoxArray &fBoxes = levels[ii]->boxArray();
        amrex::IntVect fRefRatio = levels[i]->fineRatio();

        InSituUtils::MaskCoveredCells<unsigned char>(
            pdom, cBoxes, cMap, fBoxes, fRefRatio, ng, masks[i]);
    }

    // loop over levels
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        const amrex::DistributionMapping &dMap = levels[i]->DistributionMap();

        // mask arrays for this level
        std::vector<unsigned char*> &mask = masks[i];

        // loop over boxes
        const amrex::BoxArray& ba = levels[i]->boxArray();
        unsigned int nBoxes = ba.size();

        for (unsigned int j = 0; j < nBoxes; ++j)
        {
            // skip non-local blocks
            if (dMap[j] != rank)
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
int AmrDataAdaptor::AddGhostNodesArray(vtkDataObject *mesh,
    const std::string &meshName)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::AddGhostNodesArray");

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    return 0;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName,
    int association, const std::string &arrayName)
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::AddArray");

    int rank = 0;
    MPI_Comm_rank(this->GetCommunicator(), &rank);

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

    if (!this->Internals->SimData)
    {
        SENSEI_ERROR("No simulation data")
        return -1;
    }

    if ((association != vtkDataObject::CELL) &&
        (association != vtkDataObject::POINT))
    {
        SENSEI_ERROR("Invalid association " << association)
        return -1;
    }

    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
        this->Internals->SimData->getAmrLevels();

    // find the indices of the multifab and component within for
    // the named array
    int fab = 0;
    int comp = 0;
    if (this->Internals->SimMetadata.GetIndex(arrayName, association, fab, comp))
    {
        SENSEI_ERROR("Failed to locate descriptor for "
            << sensei::VTKUtils::GetAttributesName(association)
            << " data array \"" << arrayName << "\"")
        return -1;
    }

    // loop over levels
    unsigned int nLevels = numActiveLevels(levels);
    for (unsigned int i = 0; i < nLevels; ++i)
    {
        // domain decomp
        const amrex::DistributionMapping &dmap = levels[i]->DistributionMap();

        // ghost zones
        amrex::MultiFab& state = levels[i]->get_new_data(fab);
        unsigned int ng = state.nGrow();

        if (!((association == vtkDataObject::CELL) && state.is_cell_centered()) &&
            !((association == vtkDataObject::POINT) && state.is_nodal()))
        {
            SENSEI_ERROR("association does not match MultiFAB centering")
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
        const amrex::BoxArray& ba = levels[i]->boxArray();
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
int AmrDataAdaptor::ReleaseData()
{
    sensei::TimeEvent<64> event("AmrDataAdaptor::ReleaseData");

    this->Internals->SimData = nullptr;

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
