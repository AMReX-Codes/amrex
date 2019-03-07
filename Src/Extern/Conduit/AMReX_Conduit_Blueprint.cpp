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

    // create rectilinear coordset 
    // (which also holds all implicit details needed for the topology)
    res["coordsets/coords/type"] = "rectilinear";
    res["coordsets/coords/values/x"] = DataType::float64(nx+1);
    res["coordsets/coords/values/y"] = DataType::float64(ny+1);
    
    float64_array x_coords = res["coordsets/coords/values/x"].value();
    float64_array y_coords = res["coordsets/coords/values/y"].value(); 

    float64 vx = x_min;
    for(index_t i =0; i< nx+1; i++)
    {
        x_coords[i] = vx;
        vx+=level_dx;
    }
    
    float64 vy = y_min;
    for(index_t i =0; i< ny+1; i++)
    {
        y_coords[i] = vy;
        vy+=level_dy;
    }
    
    if(dims > 2)
    {
        res["coordsets/coords/values/z"] = DataType::float64(nz+1);
        float64_array z_coords = res["coordsets/coords/values/z"].value();
        
        float64 vz = z_min;
        for(index_t i =0; i< nz+1; i++)
        {
            z_coords[i] = vz;
            vz+=level_dz;
        }
    }
    
    // create a rectilinear topology that refs our coordset
    res["topologies/topo/type"] = "rectilinear";
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
    Node &n_fields = res["fields"];
    for(int i=0; i < varnames.size(); i++)
    {
        Node &n_field = n_fields[varnames[i]];
        n_field["association"] = "element";
        n_field["topology"] = "topo";
        // const_cast is used b/c zero copy via Node::set_external 
        // requires non-const
        Real *data_ptr = const_cast<Real*>(fab.dataPtr());
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
    res.reset();
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
            patch["state/time"] = time_value;

            const FArrayBox &fab = mf[mfi];
            // create coordset and topo
            FabToBlueprintTopology(geom,fab,patch);
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
    // blueprint verify makes sure we conform to whats expected
    // for a multi-domain mesh 
    if(!blueprint::mesh::verify(res,info))
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
    else // save 1 file per mpi task + root file
    {
        //
        // For  processor case, save everything (including the 
        // blueprint index) to one file
        // 
        
        // create a folder to hold data files
        UtilCreateCleanDirectory(bp_base, true);

        std::string fname_data_base = bp_base + "/" + bp_base + "_data_";

        // save this rank's data
        std::string fname_data = amrex::Concatenate(fname_data_base,
                                                    rank,
                                                    6);
        fname_data += "." + protocol;
        relay::io::save(bp_mesh, fname_data, protocol);

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

            root["number_of_files"] = ntasks;
            root["number_of_trees"].set_int64(num_domains);
            root["file_pattern"] =  fname_data_base + "%06d." + protocol;
            root["tree_pattern"] = "domain_%06d";
            relay::io::save(root,bp_root_file,protocol);
        }

    }
}



}
