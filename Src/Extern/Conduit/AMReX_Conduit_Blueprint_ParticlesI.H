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
                        conduit::Node &res)
{
    int num_particles = ptile.GetArrayOfStructs().size();
    int struct_size   = sizeof(Particle<NStructReal, NStructInt>);

    // knowing the above, we can zero copy the x,y,z positions + id, cpu
    // and any user fields in the AOS
    
    // get the first particle's struct
    const auto &p_aos = ptile.GetArrayOfStructs()[0];

    // setup a blueprint description for the particle mesh
    
    conduit::Node &n_coords  = res["coordsets/particle_coords"];
    n_coords["type"] = "explicit";

    res["topologies/particles/coordset"] = "particle_coords";
    // create an explicit points topology
    res["topologies/particles/type"] = "unstructured";
    res["topologies/particles/elements/shape"] = "point";
    res["topologies/particles/elements/connectivity"].set(
                                    conduit::DataType::c_int(num_particles));
    int *conn = res["topologies/particles/elements/connectivity"].value();

    for(int i = 0; i < num_particles ; i++)
    {
        conn[i] = i;
    }

    //----------------------------------//
    // point locations from from aos
    //----------------------------------//

    n_coords["values/x"].set_external(const_cast<Real*>(&p_aos.m_rdata.pos[0]),
                                      num_particles,
                                      0,
                                      struct_size);
#if AMREX_SPACEDIM > 1
    n_coords["values/y"].set_external(const_cast<Real*>(&p_aos.m_rdata.pos[1]),
                                      num_particles,
                                      0,
                                      struct_size);
#endif
#if AMREX_SPACEDIM > 2
    n_coords["values/z"].set_external(const_cast<Real*>(&p_aos.m_rdata.pos[2]),
                                      num_particles,
                                      0,
                                      struct_size);
#endif

    // fields
    conduit::Node &n_fields = res["fields"];

    //----------------------------------//
    // standard integer fields from aos
    // (id, cpu)
    //----------------------------------//

    // id is the first int entry
    conduit::Node &n_f_id = n_fields["particle_id"];

    n_f_id["topology"] = "particles";
    n_f_id["association"] = "element";
    n_f_id["values"].set_external(const_cast<int*>(&p_aos.m_idata.arr[0]),
                                  num_particles,
                                  0,
                                  struct_size);

    // cpu is the second int entry
    conduit::Node &n_f_cpu = n_fields["particle_cpu"];

    n_f_cpu["topology"] = "particles";
    n_f_cpu["association"] = "element";
    n_f_cpu["values"].set_external(const_cast<int*>(&p_aos.m_idata.arr[1]),
                                   num_particles,
                                   0,
                                   struct_size);

    // -------------------------
    // user defined aos fields
    // -------------------------

    int vname_real_idx = 0;
    // struct real fields, the first set are always the particle positions
    // which we wrap above
    for (int i = AMREX_SPACEDIM; i < AMREX_SPACEDIM + NStructReal; i++)
    {
        conduit::Node &n_f = n_fields[real_comp_names[vname_real_idx]];
        n_f["topology"] = "particles";
        n_f["association"] = "element";
        n_f["values"].set_external(const_cast<Real*>(&p_aos.m_rdata.arr[i]),
                                   num_particles,
                                   0,
                                   struct_size);

        vname_real_idx++;
    }

    int vname_int_idx = 0;
    // struct int fields, first two are always id, and cpu --
    // which we wrap above
    for (int i = 2; i < 2 + NStructInt; i++)
    {
        conduit::Node &n_f = n_fields[int_comp_names[vname_int_idx]];
        n_f["topology"] = "particles";
        n_f["association"] = "element";
        n_f["values"].set_external(const_cast<int*>(&p_aos.m_idata.arr[i]),
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
        n_f["topology"] = "particles";
        n_f["association"] = "element";
        n_f["values"].set_external(const_cast<Real*>(&soa.GetRealData(i)[0]),
                                   num_particles);

        vname_real_idx++;
    }

    // array int fields
    for (int i = 0; i < NArrayInt; i++)
    {
        conduit::Node &n_f = n_fields[int_comp_names[vname_int_idx]];
        n_f["topology"] = "particles";
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
                             conduit::Node &res)
{
    res.reset();
    BL_PROFILE("ParticleContainerToBlueprint()");
    
    // validate varnames, which are used to provide field names
    // for user defined aos and soa values.

    BL_ASSERT(real_comp_names.size() == (NStructReal + NArrayReal) );
    BL_ASSERT(int_comp_names.size()  == (NStructInt + NArrayInt) );

    int num_levels = pc.maxLevel() + 1;
    int num_domains = 0;

    // loop over levels
    for (int lev = 0; lev < num_levels; lev++)
    {
        // loop over particle tiles
        for (const auto& kv : pc.GetParticles(lev)) 
        {
            // get each particle tile
            const auto& ptile = kv.second;
            
            // TODO: We need more info to provide unique domain ids 
            // across mpi tasks
            int domain_id = num_domains;

            const std::string& patch_name = amrex::Concatenate("domain_",
                                                               domain_id,
                                                               6); 
            conduit::Node &patch = res[patch_name];
            // add basic state info
            patch["state/domain_id"] = domain_id;
            ParticleTileToBlueprint(ptile,
                                    real_comp_names,
                                    int_comp_names,
                                    patch);
            num_domains++;
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
                       << info.to_json();
    }
}



}
