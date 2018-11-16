#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <sstream>

#include <AMReX_Print.H>

using namespace amrex;

namespace {

// returns coordinate in an index space with no switches
std::array<int, 4> read_node_coord (const std::string & name)
{
    int cabx, caby, cab_chas, slot, node;
    {
        std::ifstream ifs {"/proc/cray_xt/cname"};
        if (!ifs) {
            // not on a cray
            return std::array<int, 4> {0,0,0,0}; // initializer_list
        }
        char t0, t1, t2, t3, t4;
        ifs >> t0 >> cabx >> t1 >> caby >> t2 >> cab_chas >> t3 >> slot >> t4 >> node;
        assert(t0 == 'c' && t1 == '-' && t2 == 'c' && t3 == 's' && t4 == 'n');
    }

    int group = 0;
    if (name == "edison") {
        group = cabx / 2 + caby * 4; // 2 cabinets per group, 4 groups per row
        if (group > 12) { group--; } // nominal "group 12" is missing
    } else if (name == "cori-haswell") {
        group = cabx / 2;  // 2 cabinets per group
        assert(caby == 0); // only 1 row of cabinets
    } else if (name == "cori-knl") {
        group = cabx / 2 + caby * 6; // 2 cabinets per group, 6 groups per row
    } else {
        Print() << "Could not determine group!";
        std::abort();
    }
    int chas = cab_chas + 3*(cabx & 1); // 2 cabinets per group (6 chassis per group)

    return std::array<int, 4> {node, slot, chas, group};
}

std::string get_mpi_processor_name ()
{
    std::string result;
#if BL_USE_MPI
    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);
    result = std::string(name);
#endif
    return result;
}

int coord_to_id (std::array<int, 4> c)
{
    return c[0] + 4 * (c[1] + 16 * (c[2] + 6 * c[3]));
}

std::array<int, 4> id_to_coord (int id)
{
    int node = id % 4;  id /= 4;
    int slot = id % 16; id /= 16;
    int chas = id % 6;  id /= 6;
    int group = id;
    return std::array<int, 4> {node, slot, chas, group};
}

template <class T, size_t N>
std::string to_str(const std::array<T, N> & a)
{
    std::ostringstream oss;
    oss << "(";
    bool first = true;
    for (int i = 0; i < N; ++i) {
        if (!first) oss << ",";
        oss << a[i];
        first = false;
    }
    oss << ")";
    return oss.str();
}

class Machine
{
  public:
    Machine () = default;
    void init () {
        get_machine_envs();
        node_ids = get_node_ids();
    }

  private:

    const bool flag_verbose = true;

    std::string hostname;
    std::string nersc_host;
    std::string partition;
    std::string node_list;
    std::string topo_addr;

    int my_node_id;
    std::vector<int> node_ids;

    void get_machine_envs ()
    {
        hostname   = std::string(std::getenv("HOSTNAME"));
        nersc_host = std::string(std::getenv("NERSC_HOST"));
        partition  = std::string(std::getenv("SLURM_JOB_PARTITION"));
        node_list  = std::string(std::getenv("SLURM_NODELIST"));
        topo_addr  = std::string(std::getenv("SLURM_TOPOLOGY_ADDR"));

        if (flag_verbose) {
            Print() << "HOSTNAME = " << hostname << std::endl;
            Print() << "NERSC_HOST = " << nersc_host << std::endl;
            Print() << "SLURM_JOB_PARTITION = " << partition << std::endl;
            Print() << "SLURM_NODELIST = " << node_list << std::endl;
            Print() << "SLURM_TOPOLOGY_ADDR = " << topo_addr << std::endl;
        }
    }

    // get this rank's machine node ID
    int get_my_node_id ()
    {
        int result = -1;
        std::string tag = "nid";
        auto pos = topo_addr.find(tag);
        if (pos != std::string::npos) {
            result = stoi(topo_addr.substr(pos + tag.size())); // assumes format '.*nid(\d+)'
            Print() << "Got node ID from SLURM_TOPOLOGY_ADDR: " << result << std::endl;
#if BL_USE_MPI
        } else {
            auto mpi_proc_name = get_mpi_processor_name();
            Print() << "MPI_Get_processor_name: " << mpi_proc_name << std::endl;
            if (mpi_proc_name.substr(0,3) == "nid") {
                result = stoi(mpi_proc_name.substr(3)); // assumes format 'nid(\d+)'
                Print() << "Got node ID from MPI_Get_processor_name(): " << result << std::endl;
            }
#endif
        }
#ifndef NDEBUG
        auto coord = read_node_coord("cori-knl");
        int id_from_coord = coord_to_id(coord);
        AMREX_ASSERT(id_from_coord == result);
#endif
        return result;
    }

    // get all node IDs in this job, indexed by job rank
    std::vector<int> get_node_ids ()
    {
        std::vector<int> ids(ParallelDescriptor::NProcs(), 0);
#ifdef BL_USE_MPI
        int node_id = get_my_node_id();
        MPI_Allgather(&node_id, 1, ParallelDescriptor::Mpi_typemap<int>::type(),
                      ids.data(), 1, ParallelDescriptor::Mpi_typemap<int>::type(),
                      ParallelDescriptor::Communicator());
#endif
        if (flag_verbose) {
            Print() << "Rank: Node ID: Node Coord:" << std::endl;
            for (int i = 0; i < ids.size(); ++i) {
                Print() << "  " << i << ": " << ids[i] << ": " << to_str(id_to_coord(ids[i])) << std::endl;
            }
        }
        return ids;
    }

};

Machine the_machine;

}

namespace amrex {
namespace machine {

void init () { the_machine.init(); }

}}
