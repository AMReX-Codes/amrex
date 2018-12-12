#include "AMReX_AmrDataAdaptor.H"

#include "timer/Timer.h"
#include "Error.h"
#include "VTKUtils.h"

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
    std::vector<vtkDataObject*> ManagedObjects;
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
    numMeshes = 1;
    return 0;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetMeshName(unsigned int id, std::string &meshName)
{
    meshName = "mesh";
    return 0;
}

//-----------------------------------------------------------------------------
int AmrDataAdaptor::GetMesh(const std::string &meshName,
    bool structureOnly, vtkDataObject *&mesh)
{
    timer::MarkEvent("AmrDataAdaptor::GetMesh");

    mesh = nullptr;

    if (meshName != "mesh")
    {
        SENSEI_ERROR("no mesh named \"" << meshName << "\"")
        return -1;
    }

    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // get levels
    amrex::Vector<std::unique_ptr<amrex::AmrLevel>> &levels =
     this->Internals->SimData->getAmrLevels();

    unsigned int nLevels = numActiveLevels(levels);

    // initialize new vtk datasets
    vtkOverlappingAMR *amrMesh = vtkOverlappingAMR::New();
    Internals->ManagedObjects.push_back(amrMesh);
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
int AmrDataAdaptor::AddGhostNodesArray(vtkDataObject *mesh,
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
int AmrDataAdaptor::AddGhostCellsArray(vtkDataObject* mesh,
    const std::string &meshName)
{
    timer::MarkEvent("AmrDataAdaptor::AddGhostCellsArray");

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
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

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
int AmrDataAdaptor::AddArray(vtkDataObject* mesh, const std::string &meshName,
    int association, const std::string &arrayName)
{
    timer::MarkEvent("AmrDataAdaptor::AddArray");

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

//-----------------------------------------------------------------------------
int AmrDataAdaptor::ReleaseData()
{
    timer::MarkEvent("AmrDataAdaptor::ReleaseData");

    this->Internals->SimData = nullptr;

    // free up mesh objects we allocated
    size_t n = this->Internals->ManagedObjects.size();
    for (size_t i = 0; i < n; ++i)
        this->Internals->ManagedObjects[i]->Delete();
    this->Internals->ManagedObjects.clear();

    return 0;
}

}
