///////////////////////////////////////////////////////////////////////////////
//
// Template implementation of functions for Conduit Mesh Blueprint Support
// for AMReX Particle Containers
//
// This file is included in AMReX_Conduit_Blueprint.H when
// when USE_PARTICLES = TRUE
//
///////////////////////////////////////////////////////////////////////////////

#include <sstream>
#include <conduit/conduit_blueprint.hpp>
#include <conduit/conduit_relay.hpp>

namespace amrex
{
//---------------------------------------------------------------------------//
// Converts a AMReX Particle Tile into a Conduit Mesh Blueprint Hierarchy.
//---------------------------------------------------------------------------//
// Note:
// This is a helper function, it's not part of the AMReX Blueprint Interface.
//---------------------------------------------------------------------------//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleTileToBlueprint(const ParticleTile<NStructReal,
                                           NStructInt,
                                           NArrayReal,
                                           NArrayInt> &ptile,
                        const Vector<std::string> &real_comp_names,
                        const Vector<std::string> &int_comp_names,
                        conduit::Node &res,
                        const std::string &topology_name)
{
    int num_particles = ptile.GetArrayOfStructs().size();
    int struct_size   = sizeof(Particle<NStructReal, NStructInt>);

    // knowing the above, we can zero copy the x,y,z positions + id, cpu
    // and any user fields in the AOS

    // get the first particle's struct
    const auto &pstruct = ptile.GetArrayOfStructs();

    // setup a blueprint description for the particle mesh
    // create a coordinate set
    std::string coordset_name = topology_name + "_coords";
    conduit::Node &n_coords = res["coordsets"][coordset_name];
    n_coords["type"] = "explicit";

    // create an explicit points topology
    conduit::Node &n_topo = res["topologies"][topology_name];
    n_topo["coordset"] = coordset_name;
    n_topo["type"] = "unstructured";
    n_topo["elements/shape"] = "point";
    n_topo["elements/connectivity"].set(
                                    conduit::DataType::c_int(num_particles));
    int *conn = n_topo["elements/connectivity"].value();

    for(int i = 0; i < num_particles ; i++)
    {
        conn[i] = i;
    }

    //----------------------------------//
    // point locations from from aos
    //----------------------------------//

    const char* pbuf_const = reinterpret_cast<const char*>(pstruct.data());
    char* pbuf = const_cast<char*>(pbuf_const);

    ParticleReal* xp = reinterpret_cast<ParticleReal*>(pbuf); pbuf += sizeof(ParticleReal);
    n_coords["values/x"].set_external(xp,
                                      num_particles,
                                      0,
                                      struct_size);
#if AMREX_SPACEDIM > 1
    ParticleReal* yp = reinterpret_cast<ParticleReal*>(pbuf); pbuf += sizeof(ParticleReal);
    n_coords["values/y"].set_external(yp,
                                      num_particles,
                                      0,
                                      struct_size);
#endif
#if AMREX_SPACEDIM > 2
    ParticleReal* zp = reinterpret_cast<ParticleReal*>(pbuf); pbuf += sizeof(ParticleReal);
    n_coords["values/z"].set_external(zp,
                                      num_particles,
                                      0,
                                      struct_size);
#endif

    // fields
    conduit::Node &n_fields = res["fields"];

    // -----------------------------
    // user defined, real aos fields
    // -----------------------------

    int vname_real_idx = 0;
    // struct real fields, the first set are always the particle positions
    // which we wrap above
    for (int i = 0; i < NStructReal; i++)
    {
        ParticleReal* val = reinterpret_cast<ParticleReal*>(pbuf); pbuf += sizeof(ParticleReal);
        conduit::Node &n_f = n_fields[real_comp_names[vname_real_idx]];
        n_f["topology"] = topology_name;
        n_f["association"] = "element";
        n_f["values"].set_external(val,
                                   num_particles,
                                   0,
                                   struct_size);

        vname_real_idx++;
    }

    //----------------------------------//
    // standard integer fields from aos
    // (id, cpu)
    //----------------------------------//

    // id is the first int entry
    int* id = reinterpret_cast<int*>(pbuf); pbuf += sizeof(int);
    conduit::Node &n_f_id = n_fields[topology_name + "_id"];

    n_f_id["topology"] = topology_name;
    n_f_id["association"] = "element";
    n_f_id["values"].set_external(id,
                                  num_particles,
                                  0,
                                  struct_size);

    // cpu is the second int entry
    int* cpu = reinterpret_cast<int*>(pbuf); pbuf += sizeof(int);
    conduit::Node &n_f_cpu = n_fields[topology_name + "_cpu"];

    n_f_cpu["topology"] = topology_name;
    n_f_cpu["association"] = "element";
    n_f_cpu["values"].set_external(cpu,
                                   num_particles,
                                   0,
                                   struct_size);

    // --------------------------------
    // user defined, integer aos fields
    // --------------------------------

    int vname_int_idx = 0;
    for (int i = 0; i < NStructInt; i++)
    {
        int* val = reinterpret_cast<int*>(pbuf); pbuf += sizeof(int);
        conduit::Node &n_f = n_fields[int_comp_names[vname_int_idx]];
        n_f["topology"] = topology_name;
        n_f["association"] = "element";
        n_f["values"].set_external(val,
                                   num_particles,
                                   0,
                                   struct_size);
        vname_int_idx++;
    }

    // -------------------------
    // user defined soa fields
    // -------------------------

     const auto &soa = ptile.GetStructOfArrays();

    // for soa entries, we can use standard strides,
    // since these are contiguous arrays

    // array real fields
    for (int i = 0; i < NArrayReal; i++)
    {
        conduit::Node &n_f = n_fields[real_comp_names[vname_real_idx]];
        n_f["topology"] = topology_name;
        n_f["association"] = "element";
        n_f["values"].set_external(const_cast<ParticleReal*>(&soa.GetRealData(i)[0]),
                                   num_particles);

        vname_real_idx++;
    }

    // array int fields
    for (int i = 0; i < NArrayInt; i++)
    {
        conduit::Node &n_f = n_fields[int_comp_names[vname_int_idx]];
        n_f["topology"] = topology_name;
        n_f["association"] = "element";
        n_f["values"].set_external(const_cast<int*>(&soa.GetIntData(i)[0]),
                                   num_particles);

        vname_int_idx++;
    }
}

//---------------------------------------------------------------------------//
// Converts a AMReX Particle Container into a Conduit Mesh Blueprint Hierarchy.
//---------------------------------------------------------------------------//
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void 
ParticleContainerToBlueprint(const ParticleContainer<NStructReal,
                                                     NStructInt,
                                                     NArrayReal,
                                                     NArrayInt> &pc,
                             const Vector<std::string> &real_comp_names,
                             const Vector<std::string> &int_comp_names,
                             conduit::Node &res,
                             const std::string &topology_name)
{
    BL_PROFILE("ParticleContainerToBlueprint()");
    
    // validate varnames, which are used to provide field names
    // for user defined aos and soa values.

    BL_ASSERT(real_comp_names.size() == (NStructReal + NArrayReal) );
    BL_ASSERT(int_comp_names.size()  == (NStructInt + NArrayInt) );

    int num_levels = pc.maxLevel() + 1;
    int num_domains = 0;

    // loop over levels

    int rank   = ParallelDescriptor::MyProc();
    int nprocs = ParallelDescriptor::NProcs();
    
    using MyParConstIter = ParConstIter<NStructReal, NStructInt, NArrayReal, NArrayInt>;

    //
    // blueprint expects unique ids for each domain published
    // for particle data, this means we need a unique id for each
    // tile across levels.
    //
    // to generate unique domain ids, we use a counting loop to
    // calculate a per-level, per-rank offset
    //
    
    // this holds our unique offset for each level for the current rank
    std::vector<int> my_lvl_offsets;
    my_lvl_offsets.resize(num_levels);

    // these hold the values we need to do a scan to calc offsets
    std::vector<int> within_lvl_counts;
    std::vector<int> within_lvl_offsets;
    within_lvl_counts.resize(nprocs);
    within_lvl_offsets.resize(nprocs);

    int total_num_domains = 0;
    for (int lev = 0; lev < num_levels; ++lev) 
    {
        // count how many tiles this mpi task has
        long level_num_domains = 0;
        for (MyParConstIter pti(pc, lev); pti.isValid(); ++pti)
        {
            level_num_domains+=1;
        }

        // clear the within_lvl_counts and offsets data
        for(int rank_idx = 0; rank_idx < nprocs ; rank_idx++)
        {
            within_lvl_counts[rank_idx]  = 0;
            within_lvl_offsets[rank_idx] = 0;
        }
        
        // set number of domains for current rank at current level
        within_lvl_counts[rank] = level_num_domains;
        // mpi sum reduce to get all rank's counts for this level
        ParallelDescriptor::ReduceIntSum(&within_lvl_counts[0],
                                         nprocs);
            
        long level_total_num_domains = 0;
        // scan to calc offset and total number of domains
        for(int rank_idx = 0; rank_idx < nprocs ; rank_idx++)
        {
            within_lvl_offsets[rank_idx] = level_total_num_domains;
            level_total_num_domains += within_lvl_counts[rank_idx];
        }

        // the offset for this rank at this level is the current 
        // total number of domains + the offset for this level
        my_lvl_offsets[lev] = within_lvl_offsets[rank] + total_num_domains;

        total_num_domains += level_num_domains;
    }

    for (int lev = 0; lev < num_levels; ++lev) 
    {
        // get our unique global offset for this rank at this level
        int lvl_offset = my_lvl_offsets[lev];
        // keep a counter of local # of tiles
        int num_lvl_local_tiles = 0;
        for (MyParConstIter pti(pc, lev); pti.isValid(); ++pti)
        {
            // this will give us a unique id across levels for each tile
            int domain_id = lvl_offset + num_lvl_local_tiles;

            const auto& ptile = pti.GetParticleTile();
            const std::string& patch_name = amrex::Concatenate("domain_",
                                                               domain_id,
                                                               6);
            conduit::Node &patch = res[patch_name];
            // add basic state info
            patch["state/domain_id"] = domain_id;
            ParticleTileToBlueprint(ptile,
                                    real_comp_names,
                                    int_comp_names,
                                    patch,
                                    topology_name);
            num_lvl_local_tiles += 1;
        }
    }

    conduit::Node info;
    // blueprint verify makes sure we conform to whats expected
    // for a multi-domain mesh 
    if(!conduit::blueprint::mesh::verify(res,info))
    {
        // ERROR -- doesn't conform to the mesh blueprint
        // show what went wrong
        amrex::Print() << "ERROR: Conduit Mesh Blueprint Verify Failed!\n"
                       << info.to_yaml();
    }
}



}
