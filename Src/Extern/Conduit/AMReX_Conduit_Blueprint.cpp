///////////////////////////////////////////////////////////////////////////////
//
// Conduit Mesh Blueprint Support for AMReX Meshes
//
// This file is only compiled when USE_CONDUIT = TRUE
//
///////////////////////////////////////////////////////////////////////////////

#include <AMReX_Conduit_Blueprint.H>

#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>
using namespace conduit;

namespace amrex {

//---------------------------------------------------------------------------//
// Creates the Blueprint nesting relationships
//---------------------------------------------------------------------------//
bool Nestsets(const int level,
              const int n_levels,
              const FArrayBox &fab,
              const Vector<const BoxArray*> box_arrays,
              const Vector<IntVect> &ref_ratio,
              const Vector<int> &domain_offsets,
              conduit::Node &nestset)
{
    nestset.reset();
    nestset["association"] = "element";
    nestset["topology"] = "topo";

    const int dims = BL_SPACEDIM;
    const Box &box = fab.box();

    bool valid = false;
    if(level > 0)
    {
      // check for parents
      std::vector<std::pair<int,Box> > isects
        = box_arrays[level-1]->intersections(amrex::coarsen(box, ref_ratio[level-1]));

      for(int b = 0; b < isects.size(); ++b)
      {
        valid = true;
        // get parent box in terms of this level
        Box parent = amrex::refine(isects[b].second, ref_ratio[level-1]);
        Box overlap = box & parent;
        int parent_id = isects[b].first + domain_offsets[level-1];

        const std::string& w_name = amrex::Concatenate("window_",
                                                        parent_id,
                                                        4);
        conduit::Node &window = nestset["windows/"+w_name];
        window["domain_id"] = parent_id;
        window["domain_type"] = "parent";
        // box coordinates are global to the level,
        // but the the window is local to this box so
        // subtract the current box origin
        window["origin/i"] = overlap.smallEnd()[0] - box.smallEnd()[0];
        window["origin/j"] = overlap.smallEnd()[1] - box.smallEnd()[1];
        if(dims == 3)
        {
            window["origin/k"] = overlap.smallEnd()[2] - box.smallEnd()[2];
        }
        window["dims/i"] = overlap.size()[0];
        window["dims/j"] = overlap.size()[1];
        if(dims == 3)
        {
            window["dims/k"] = overlap.size()[2];
        }
        window["ratio/i"] = ref_ratio[level-1][0];
        window["ratio/j"] = ref_ratio[level-1][1];
        if(dims == 3)
        {
            window["ratio/k"] = ref_ratio[level-1][2];
        }
      }
    }

    if(level < n_levels - 1)
    {
      // check for children
      std::vector<std::pair<int,Box> > isects
        = box_arrays[level+1]->intersections(amrex::refine(box, ref_ratio[level]));

      for(int b = 0; b < isects.size(); ++b)
      {
        valid = true;
        // get child box in terms of this level
        Box child = amrex::coarsen(isects[b].second, ref_ratio[level]);
        int child_id = isects[b].first + domain_offsets[level+1];
        Box overlap = box & child;

        const std::string& w_name = amrex::Concatenate("window_",
                                                        child_id,
                                                        4);

        conduit::Node &window = nestset["windows/"+w_name];
        window["domain_id"] = child_id;
        window["domain_type"] = "child";
        // box coordinates are global to the level,
        // but the the window is local to this box so
        // subtract the current box origin
        window["origin/i"] = overlap.smallEnd()[0] - box.smallEnd()[0];
        window["origin/j"] = overlap.smallEnd()[1] - box.smallEnd()[1];
        if(dims == 3)
        {
            window["origin/k"] = overlap.smallEnd()[2] - box.smallEnd()[2];
        }
        window["dims/i"] = overlap.size()[0];
        window["dims/j"] = overlap.size()[1];
        if(dims == 3)
        {
            window["dims/k"] = overlap.size()[2];
        }
        window["ratio/i"] = ref_ratio[level+1][0];
        window["ratio/j"] = ref_ratio[level+1][1];
        if(dims == 3)
        {
            window["ratio/k"] = ref_ratio[level+1][2];
        }
      }
    }
    return valid;
}

//---------------------------------------------------------------------------//
// Creates a Rectilinear Blueprint Topology from the Geom and Fab.
//---------------------------------------------------------------------------//
void FabToBlueprintTopology(const Geometry& geom,
                            const FArrayBox& fab,
                            Node &res)
{
    int dims = BL_SPACEDIM;

    // get the details of the entire level from geom
    amrex::Box level_box = geom.Domain();

    int level_nx = level_box.size()[0];
    int level_ny = level_box.size()[1];
    int level_nz = dims > 2 ? level_box.size()[2] : 0;

    float64 level_x_min = geom.ProbLo(0);
    float64 level_x_max = geom.ProbHi(0);

    float64 level_y_min = geom.ProbLo(1);
    float64 level_y_max = geom.ProbHi(1);

    float64 level_z_min = dims > 2 ? geom.ProbLo(2) : 0;
    float64 level_z_max = dims > 2 ? geom.ProbHi(2) : 0;

    // geom.CellSize()[i] == (level_x_max - level_x_min) / float64(level_nx);
    float64 level_dx = geom.CellSize()[0];

    // geom.CellSize()[j] == (level_y_max - level_y_min) / float64(level_ny);
    float64 level_dy = geom.CellSize()[1];

    // geom.CellSize()[k] == (level_z_max - level_z_min) / float64(level_nz)
    float64 level_dz = dims > 2 ? geom.CellSize()[2] : 0.0;

    // now extract the FAB details
    const amrex::Box &fab_box = fab.box();


    int i_min = fab_box.smallEnd(0);
    int i_max = fab_box.bigEnd(0);

    int j_min = fab_box.smallEnd(1);
    int j_max = fab_box.bigEnd(1);

    int k_min = dims > 2 ? fab_box.smallEnd(2) : 0;
    int k_max = dims > 2 ? fab_box.bigEnd(2) : 0;


    int nx = (i_max - i_min + 1);
    int ny = (j_max - j_min + 1);
    int nz = dims > 2 ? (k_max - k_min +1) : 1;

    float64 x_min = level_x_min + level_dx * i_min;
    float64 x_max = level_x_min + level_dx * i_max;

    float64 y_min = level_y_min + level_dy * j_min;
    float64 y_max = level_y_min + level_dy * j_max;

    float64 z_min = dims > 2 ? level_z_min + level_dz * k_min : 0.0;
    float64 z_max = dims > 2 ? level_z_min + level_dz * k_max : 0.0;

    // create uniform coordset
    // (which also holds all implicit details needed for the topology)
    res["coordsets/coords/type"] = "uniform";
    res["coordsets/coords/dims/i"] = nx+1;
    res["coordsets/coords/dims/j"] = ny+1;

    res["coordsets/coords/spacing/dx"] = level_dx;
    res["coordsets/coords/spacing/dy"] = level_dy;

    res["coordsets/coords/origin/x"] = x_min;
    res["coordsets/coords/origin/y"] = y_min;

    if(dims > 2)
    {
      res["coordsets/coords/dims/k"] = nz+1;
      res["coordsets/coords/spacing/dz"] = level_dz;
      res["coordsets/coords/origin/z"] = z_min;
    }

    // create a rectilinear topology that refs our coordset
    res["topologies/topo/type"] = "uniform";
    res["topologies/topo/coordset"] = "coords";

    // add logical elements origin
    res["topologies/topo/elements/origin/i0"] = i_min;
    res["topologies/topo/elements/origin/j0"] = j_min;
    if( dims > 2)
    {
        res["topologies/topo/elements/origin/k0"] = k_min;
    }

}

//---------------------------------------------------------------------------//
// Creates a Blueprint Field that identifies the cells in the ghost or "grow"
// region of the fab.
//
// The values indicate:
//    0 = normal cell
//    1 = ghost cell
//
//---------------------------------------------------------------------------//
void AddFabGhostIndicatorField (const FArrayBox& fab,
                                int ngrow,
                                Node &res)
{
    Node &n_field = res["fields/ghost_indicator"];
    n_field["association"] = "element";
    n_field["topology"] = "topo";
    n_field["values"].set(DataType::float64(fab.box().numPts()));
    float64_array vals_array = n_field["values"].value();

    int dims = BL_SPACEDIM;

    // find the FAB details
    const amrex::Box &fab_box = fab.box();

    int i_min = fab_box.smallEnd(0);
    int i_max = fab_box.bigEnd(0);

    int j_min = fab_box.smallEnd(1);
    int j_max = fab_box.bigEnd(1);

    int k_min = dims > 2 ? fab_box.smallEnd(2) : 0;
    int k_max = dims > 2 ? fab_box.bigEnd(2) : 0;

    int ni = (i_max - i_min + 1);
    int nj = (j_max - j_min + 1);
    int nk = dims > 2 ? (k_max - k_min +1) : 1;

    // ngrow is the number of ghosts on each side of the box
    // set indicator to 1 if the zone is a ghost
    int idx = 0;
    for(int k=0; k < nk; k++)
    {
        for(int j=0; j < nj; j++)
        {
            for(int i=0; i < ni; i++)
            {
                if(  (i < ngrow ) || ( (ni - ngrow -1) < i ) ||
                     (j < ngrow ) || ( (nj - ngrow -1) < j ) ||
                     ( (dims >2 )  &&
                       ( (k < ngrow ) || ( (nk - ngrow -1) < k ) )
                     )
                  )
                {
                    vals_array[idx] = 1;
                }
                idx++;
            }
        }
    }
}

//---------------------------------------------------------------------------//
// Wraps the fab variables into Blueprint Fields.
//---------------------------------------------------------------------------//
void FabToBlueprintFields (const FArrayBox& fab,
                           const Vector<std::string>& varnames,
                           Node &res)
{
    // make sure we are not asking for more components than exist.
    BL_ASSERT(varnames.size() <= fab->nComp());

    Node &n_fields = res["fields"];
    for(int i=0; i < varnames.size(); i++)
    {
        Node &n_field = n_fields[varnames[i]];
        n_field["association"] = "element";
        n_field["topology"] = "topo";
        //
        // const_cast is used b/c zero copy via Node::set_external
        // requires non-const
        //
        // the field data values are the i'th component
        //
        Real *data_ptr = const_cast<Real*>(fab.dataPtr(i));
        n_field["values"].set_external(data_ptr,fab.box().numPts());
    }
}

//---------------------------------------------------------------------------//
// Converts a single level AMReX mesh into a Conduit Mesh Blueprint Hierarchy.
//---------------------------------------------------------------------------//
void
SingleLevelToBlueprint (const MultiFab& mf,
                        const Vector<std::string>& varnames,
                        const Geometry& geom,
                        Real time_value,
                        int level_step,
                        Node &res)
{
    // prepare inputs  so we can call the multi-level case
    Vector<const MultiFab*> mfarr(1,&mf);
    Vector<Geometry> geomarr(1,geom);
    Vector<int> level_steps(1,level_step);
    Vector<IntVect> ref_ratio;

    // call multi-level case
    MultiLevelToBlueprint(1,
                          mfarr,
                          varnames,
                          geomarr,
                          time_value,
                          level_steps,
                          ref_ratio,
                          res);
}

//---------------------------------------------------------------------------//
// Converts a AMReX AMR mesh into a Conduit Mesh Blueprint Hierarchy.
//---------------------------------------------------------------------------//
void
MultiLevelToBlueprint (int n_levels,
                       const Vector<const MultiFab*>& mfs,
                       const Vector<std::string>& varnames,
                       const Vector<Geometry>& geoms,
                       Real time_value,
                       const Vector<int>& level_steps,
                       const Vector<IntVect>& ref_ratio,
                       Node &res)
{
    BL_PROFILE("MultiLevelToBlueprint()");

    BL_ASSERT(n_levels <= mfs.size());
    BL_ASSERT(n_levels <= geoms.size());
    BL_ASSERT(n_levels <= ref_ratio.size()+1);
    BL_ASSERT(n_levels <= level_steps.size());
    BL_ASSERT(mfs[0]->nComp() == varnames.size());

    int finest_level = n_levels-1;

    int num_levels = geoms.size();

    // get mpi rank and # of tasks
    int rank   = ParallelDescriptor::MyProc();
    int ntasks = ParallelDescriptor::NProcs();

    Vector<const BoxArray*> box_arrays;
    Vector<int> box_offsets;
    for(int i = 0; i < num_levels; i++)
    {
      const BoxArray &boxs = mfs[i]->boxArray();
      box_arrays.push_back(&boxs);
      if(i == 0)
      {
        box_offsets.push_back(0);
      }
      else
      {
        box_offsets[i] = box_offsets[i-1] + mfs[i]->size();
      }
    }

    int num_domains = 0;
    for(int i = 0; i < num_levels; i++)
    {
        //
        // Geometry represents the physical and logical space of an entire level.
        //
        // Multifab contains the patches or blocks for this level.
        // In Blueprint speak, Each Multifab contains several domains and
        // each fab has "components" which map to a Blueprint field.

        const Geometry &geom = geoms[i];
        const MultiFab &mf = *mfs[i];

        // ngrow tells us how many layers of ghosts
        int ngrow = mf.nGrow();

        // mfiter allows us to iterate over local patches
        for(MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            // domain_id is mfi.index + all patches on lower levels
            int domain_id = mfi.index() + num_domains;
            const std::string& patch_name = amrex::Concatenate("domain_",
                                                               domain_id,
                                                               6);

            Node &patch = res[patch_name];
            // add basic state info
            patch["state/domain_id"] = domain_id;
            patch["state/cycle"] = level_steps[0];
            patch["state/time"] = time_value;

            const FArrayBox &fab = mf[mfi];

            // create coordset and topo
            FabToBlueprintTopology(geom,fab,patch);
            // add the nesting relationship
            if(num_levels > 1)
            {
                conduit::Node nestset;
                bool valid = Nestsets(i, n_levels, fab, box_arrays, ref_ratio, box_offsets, nestset);
                if(valid)
                {
                    patch["nestsets/nest"].set(nestset);
                }
            }
            // add fields
            FabToBlueprintFields(fab,varnames,patch);

            // add ghost indicator if the fab has ghost cells
            if(ngrow > 0)
            {
                AddFabGhostIndicatorField(fab,ngrow,patch);
            }
        }
        num_domains += mf.size();
    }

    Node info;
    // if we have mesh data, use blueprint verify
    // to make sure we conform to whats expected
    // for a multi-domain mesh

    if(!res.dtype().is_empty() &&
       !blueprint::mesh::verify(res,info))
    {
        // ERROR -- doesn't conform to the mesh blueprint
        // show what went wrong
        amrex::Print() << "ERROR: Conduit Mesh Blueprint Verify Failed!\n"
                       << info.to_json();
    }

}

//---------------------------------------------------------------------------//
// Write a Conduit Mesh Blueprint Hierarchy to a set of files that can
// be viewed in Visit using the Blueprint plugin.
//---------------------------------------------------------------------------//
void WriteBlueprintFiles (const conduit::Node &bp_mesh,
                          const std::string &fname_base,
                          int   step,
                          const std::string &protocol)
{
    BL_PROFILE("WriteBlueprintData()");

    // find the # of domains
    long num_domains = (long)bp_mesh.number_of_children();
    ParallelDescriptor::ReduceLongSum(num_domains);

    // get numer of mpi tasks and this task's rank
    int rank   = ParallelDescriptor::MyProc();
    int ntasks = ParallelDescriptor::NProcs();

    // gen base file name that includes padded step #
    const std::string& bp_base = amrex::Concatenate(fname_base,
                                                    step,
                                                    5);

    std::string bp_root_file = bp_base + ".blueprint_root";

    //
    // For the 1 processor case, save everything (including the
    // blueprint index) to one file
    //
    if(ntasks == 1 )  // save everything to one file
    {
        // we don't want to modify the input tree,
        // so we  zero-copy the mesh blueprint tree to
        // and add the index info for i/o
        Node root;
        root.set_external(bp_mesh);
        // gen index from first domain
        blueprint::mesh::generate_index(bp_mesh[0],
                                    "",
                                    num_domains,
                                    root["blueprint_index/mesh"]);
        // add final root information
        root["protocol/name"] = protocol;
        root["protocol/version"] = "0.4.0";

        root["number_of_files"] = ntasks; // 1
        root["number_of_trees"].set_int64(num_domains);
        root["file_pattern"] = bp_root_file;
        root["tree_pattern"] = "domain_%06d";
        relay::io::save(root,bp_root_file,protocol);
    }
    else // save 1 per file domain + root file
    {
        // one file per domain, until blueprint clients
        // get generalized

        // create a folder to hold data files
        UtilCreateCleanDirectory(bp_base, true);

        std::string fname_data_base = bp_base + "/" + bp_base + "_data_";

        // save each domain to its own file
        NodeConstIterator itr = bp_mesh.children();
        while(itr.has_next())
        {
            const Node &n = itr.next();

            // fetch the domain id
            int64 domain_id = n["state/domain_id"].to_int64();

            // generate file for this domain
            std::string fname_data = amrex::Concatenate(fname_data_base,
                                                        domain_id,
                                                        6);
            fname_data += "." + protocol;

            std::string domain_path = amrex::Concatenate(":domain_",
                                                         domain_id,
                                                         6);
            // save this domain's data
            relay::io::save(n,
                            fname_data + domain_path,
                            protocol);
        }

        // generate root file on the MPI Task designated as
        // the AMReX IO Processor
        if(rank == ParallelDescriptor::IOProcessorNumber())
        {
            // gen index from first local domain
            Node root;
            blueprint::mesh::generate_index(bp_mesh[0],
                                            "",
                                            num_domains,
                                            root["blueprint_index/mesh"]);

            // add final root information
            root["protocol/name"] = protocol;
            root["protocol/version"] = "0.4.0";

            root["number_of_files"] = num_domains;
            root["number_of_trees"].set_int64(num_domains);
            root["file_pattern"] =  fname_data_base + "%06d." + protocol;
            root["tree_pattern"] = "domain_%06d";
            relay::io::save(root,bp_root_file,protocol);
        }
    }
}


}
