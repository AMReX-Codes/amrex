#include <string>
#include <vector>
#include <array>
#include <fstream>
#include <cassert>
#include <sstream>
#include <map>
#include <unordered_map>

#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>
#include <AMReX_ParallelReduce.H>
#include <AMReX_Utility.H>
#include <AMReX_Machine.H>

using namespace amrex;

namespace {

struct DoubleInt {
    double d;
    int i;
};

using Coord = Array<int, 4>;

// returns coordinate in an index space with no switches
// for dragonfly network
Coord read_df_node_coord (const std::string & name)
{
    int cabx, caby, cab_chas, slot, node;
    {
        std::ifstream ifs {"/proc/cray_xt/cname"};
        if (!ifs) {
            // not on a cray
            return Coord {0,0,0,0}; // initializer_list
        }
        char t0, t1, t2, t3, t4;
        ifs >> t0 >> cabx >> t1 >> caby >> t2 >> cab_chas >> t3 >> slot >> t4 >> node;
        assert(t0 == 'c' && t1 == '-' && t2 == 'c' && t3 == 's' && t4 == 'n');
    }

    int group = 0;
    if (name == "edison") {
        group = cabx / 2 + caby * 4; // 2 cabinets per group, 4 groups per row
        if (group > 12) { group--; } // nominal "group 12" is missing
    } else if (name == "cori") {
        group = cabx / 2 + caby * 6; // 2 cabinets per group, 6 groups per row
    } else {
        Print() << "Could not determine group!";
        std::abort();
    }
    int chas = cab_chas + 3*(cabx & 1); // 2 cabinets per group (6 chassis per group)

    return Coord {node, slot, chas, group};
}

std::string get_mpi_processor_name ()
{
    std::string result;
#ifdef BL_USE_MPI
    int len;
    char name[MPI_MAX_PROCESSOR_NAME];
    MPI_Get_processor_name(name, &len);
    result = std::string(name);
#endif
    return result;
}

// assumes groups are in 4x16x6 configuration
int df_coord_to_id (const Coord & c)
{
    return c[0] + 4 * (c[1] + 16 * (c[2] + 6 * c[3]));
}

// assumes groups are in 4x16x6 configuration
Coord df_id_to_coord (int id)
{
    int node = id % 4;  id /= 4;
    int slot = id % 16; id /= 16;
    int chas = id % 6;  id /= 6;
    int group = id;
    return Coord {node, slot, chas, group};
}

template <class T, size_t N>
std::string to_str(const Array<T, N> & a)
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

template <class T>
std::string to_str(const Vector<T> & v)
{
    std::ostringstream oss;
    oss << "(";
    bool first = true;
    for (int i = 0; i < v.size(); ++i) {
        if (!first) oss << ",";
        oss << v[i];
        first = false;
    }
    oss << ")";
    return oss.str();
}

Vector<int> get_subgroup_ranks ()
{
    int rank_n = ParallelContext::NProcsSub();
    Vector<int> lranks(rank_n);
    for (int i = 0; i < rank_n; ++i) {
        lranks[i] = i;
    }

    Vector<int> granks(rank_n);
    ParallelContext::local_to_global_rank(granks.data(), lranks.data(), rank_n);
    return granks;
}

int pair_n (int x) {
    return x*(x-1)/2;
}

int df_dist (const Coord & a, const Coord & b)
{
    if (a[3] != b[3]) {
        // large penalty for traversing across groups
        return 20;
    } else {
        // same group
        int slot_diff = (a[1] != b[1] ? 1 : 0);
        int chas_diff = (a[2] != b[2] ? 1 : 0);
        if (slot_diff + chas_diff == 0 && a[0] == b[0]) {
            // same node
            return 0;
        } else {
            // add 2 for first and last node-to-switch hops
            return 2 + slot_diff + chas_diff;
        }
    }
}

Coord id_to_coord (int id)
{
    // TODO: implement support for other types of networks
    return df_id_to_coord(id);
}

int dist (const Coord & a, const Coord & b)
{
    // TODO: implement support for other types of networks
    return df_dist(a, b);
}

struct Candidate
{
    int id;
    Coord coord;
    // how many ranks on this node
    int rank_n = 0;
    // sum of pairwise rank distances from the candidate node to already chosen nodes
    int sum_dist = 0;

    Candidate () = default;
    Candidate (int i) : id(i), coord(id_to_coord(id)) {}
};

class NeighborhoodCache
{
  public:
    void add (uint64_t key, Vector<int> val) {
        AMREX_ASSERT(cache.count(key) == 0);
        cache[key] = std::move(val);
    }
    bool get (uint64_t key, Vector<int> & val) {
        bool result = cache.count(key) > 0;
        if (result) {
            val = cache.at(key);
        }
        return result;
    }

    // result is dependent on both the current set of ranks
    // and the size of the neighborhood desired
    uint64_t hash (const Vector<int> & cur_ranks, int nbh_rank_n) {
        auto result = hash_vector(cur_ranks);
        hash_combine(result, nbh_rank_n);
        return result;
    }

  private:
    std::unordered_map<uint64_t, Vector<int>> cache;
};

class Machine
{
  public:
    Machine () {
        get_params();
        get_machine_envs();
        node_ids = get_node_ids();
    }

    // find a compact neighborhood of size rank_n in the current ParallelContext subgroup
    Vector<int> find_best_nbh (int nbh_rank_n, bool flag_local_ranks)
    {
#ifdef BL_USE_MPI
        BL_PROFILE("Machine::find_best_nbh()");

        auto sg_g_ranks = get_subgroup_ranks();
        auto sg_rank_n = sg_g_ranks.size();
        if (flag_verbose) {
            Print() << "Machine::find_best_nbh(): called for " << nbh_rank_n
                    << " of " << sg_rank_n << " ranks" << std::endl;
        }

        Vector<int> result;
        auto key = nbh_cache.hash(sg_g_ranks, nbh_rank_n);
        if (nbh_cache.get(key, result)) {
            if (flag_verbose) {
                Print() << "Machine::find_best_nbh(): found neighborhood in cache" << std::endl;
            }
        } else {
            // get node IDs of current subgroup
            Vector<int> sg_node_ids(sg_rank_n);
            std::unordered_map<int, std::vector<int>> node_ranks;
            for (int i = 0; i < sg_rank_n; ++i) {
                AMREX_ASSERT(sg_g_ranks[i] >= 0 && sg_g_ranks[i] < node_ids.size());
                sg_node_ids[i] = node_ids[sg_g_ranks[i]];
                if (flag_local_ranks) {
                    node_ranks[sg_node_ids[i]].push_back(i);
                } else {
                    node_ranks[sg_node_ids[i]].push_back(sg_g_ranks[i]);
                }
            }

            if (flag_very_verbose) {
                Print() << "SubRank: GloRank: Node ID: Node Coord:" << std::endl;
                for (int i = 0; i < sg_rank_n; ++i) {
                    Print() << "  " << i << ": " << sg_g_ranks[i] << ": " << sg_node_ids[i]
                            << ": " << to_str(id_to_coord(sg_node_ids[i])) << std::endl;
                }
            }

            Vector<int> local_nbh;
            double score;
            auto rank_me = ParallelContext::MyProcSub();
            tie(local_nbh, score) = search_local_nbh(rank_me, sg_node_ids, nbh_rank_n);

            if (flag_verbose) {
                Vector<int> base_nbh;
                double base_score;
                tie(base_nbh, base_score) = baseline_score(sg_node_ids, nbh_rank_n);

                Print() << "Baseline neighborhood: " << to_str(base_nbh) << ", score = " << base_score << std::endl;
                Print() << "Rank 0's neighborhood: " << to_str(local_nbh) << ", score = " << score << std::endl;
            }

            // determine the best neighborhood among ranks
            DoubleInt my_score_with_id {score, rank_me}, min_score_with_id;
            MPI_Allreduce(&my_score_with_id, &min_score_with_id, 1, MPI_DOUBLE_INT, MPI_MINLOC, ParallelContext::CommunicatorSub());
            double winner_score = min_score_with_id.d;
            int    winner_rank  = min_score_with_id.i;

            // broadcast the best hood from winner rank to everyone
            int local_nbh_size = local_nbh.size();
            MPI_Bcast(&local_nbh_size, 1, MPI_INT, winner_rank, ParallelContext::CommunicatorSub());
            local_nbh.resize(local_nbh_size);
            MPI_Bcast(local_nbh.data(), local_nbh.size(), MPI_INT, winner_rank, ParallelContext::CommunicatorSub()); 

            std::sort(local_nbh.begin(), local_nbh.end());
            if (flag_verbose) {
                Print() << "Winning neighborhood: " << winner_rank << ": " << to_str(local_nbh)
                        << ", score = " << winner_score << std::endl;
            }

            result.reserve(nbh_rank_n);
            for (int i = 0; i < local_nbh.size(); ++i) {
                for (auto rank : node_ranks.at(local_nbh[i])) {
                    if (result.size() < nbh_rank_n) {
                        result.push_back(rank);
                    }
                }
            }
            nbh_cache.add(key, result);
        }

        if (flag_very_verbose) {
            Print() << "Ranks in neighborhood: " << to_str(result) << std::endl;
        }

        return result;
#else
        return Vector<int>(nbh_rank_n, 0);
#endif
    }

  private:

    std::string hostname;
    std::string nersc_host;
    std::string partition;
    std::string node_list;
    std::string topo_addr;

    int flag_verbose = 0;
    int flag_very_verbose = 0;
    bool flag_nersc_df;
    int my_node_id;
    Vector<int> node_ids;

    NeighborhoodCache nbh_cache;

    void get_params ()
    {
        ParmParse pp("machine");
        pp.query("verbose", flag_verbose);
        pp.query("very_verbose", flag_very_verbose);
    }

    std::string get_env_str (std::string env_key)
    {
        std::string result;
        auto val_c_str = std::getenv(env_key.c_str());
        if (val_c_str) {
            result = std::string(val_c_str);
        }
        return result;
    }

    void get_machine_envs ()
    {
        hostname   = get_env_str("HOSTNAME");
        nersc_host = get_env_str("NERSC_HOST");
#ifdef AMREX_USE_CUDA
        flag_nersc_df = false;
#else
        flag_nersc_df = (nersc_host == "edison" ||
                         nersc_host == "cori" ||
                         nersc_host == "saul");
#endif

        if (flag_nersc_df) {
            partition  = get_env_str("SLURM_JOB_PARTITION");
            node_list  = get_env_str("SLURM_NODELIST");
            topo_addr  = get_env_str("SLURM_TOPOLOGY_ADDR");

            if (flag_verbose) {
                Print() << "HOSTNAME = " << hostname << std::endl;
                Print() << "NERSC_HOST = " << nersc_host << std::endl;
                Print() << "SLURM_JOB_PARTITION = " << partition << std::endl;
                Print() << "SLURM_NODELIST = " << node_list << std::endl;
                Print() << "SLURM_TOPOLOGY_ADDR = " << topo_addr << std::endl;
            }
        }
    }

    // get this rank's machine node ID
    int get_my_node_id ()
    {
        int result = -1;
        if (flag_nersc_df) {
            std::string tag = "nid";
            auto pos = topo_addr.find(tag);
            if (pos != std::string::npos) {
                result = stoi(topo_addr.substr(pos + tag.size())); // assumes format ".*nid(\d+)"
                if (flag_verbose) {
                    Print() << "Got node ID from SLURM_TOPOLOGY_ADDR: " << result << std::endl;
                }
#ifdef BL_USE_MPI
            } else {
                auto mpi_proc_name = get_mpi_processor_name();
                Print() << "MPI_Get_processor_name: " << mpi_proc_name << std::endl;
                pos = mpi_proc_name.find(tag);
                if (pos != std::string::npos) {
                    result = stoi(mpi_proc_name.substr(pos + tag.size())); // assumes format ".*nid(\d+)"
                    if (flag_verbose) {
                        Print() << "Got node ID from MPI_Get_processor_name(): " << result << std::endl;
                    }
                }
#endif
            }

            // check result
            AMREX_ALWAYS_ASSERT(result != -1);
#ifndef NDEBUG
            auto coord = read_df_node_coord(nersc_host);
            int id_from_coord = df_coord_to_id(coord);
            AMREX_ASSERT(id_from_coord == result);
#endif
        } else {
            result = 0;
        }

        return result;
    }

    // get all node IDs in this job, indexed by job rank
    // this is collective over ALL ranks in the job
    Vector<int> get_node_ids ()
    {
        int node_id = -1;
        Vector<int> ids(ParallelDescriptor::NProcs(), 0);
#ifdef BL_USE_MPI
        node_id = get_my_node_id();
        ParallelAllGather::AllGather(node_id, ids.data(), ParallelContext::CommunicatorAll());
#endif
        if (flag_verbose) {
            std::map<int, Vector<int>> node_ranks;
            for (int i = 0; i < ids.size(); ++i) {
                node_ranks[ids[i]].push_back(i);
            }
            Print() << "Node ID: Node Coord: Ranks:" << std::endl;
            for (const auto & p : node_ranks) {
                Print() << "  " << p.first << ": " << to_str(id_to_coord(p.first))
                        << ": " << to_str(p.second) << std::endl;
            }
        }
        return ids;
    }

    // do a local search starting at current node
    std::pair<Vector<int>, double>
    baseline_score(const Vector<int> & sg_node_ids, int nbh_rank_n)
    {
        AMREX_ASSERT(sg_node_ids.size() > 0 && nbh_rank_n > 0 &&
                     nbh_rank_n <= sg_node_ids.size());

        // construct map of node candidates to select
        std::map<int, Candidate> cand_map;
        for (int i = 0; i < nbh_rank_n; ++i) {
            auto node_id = sg_node_ids[i];
            if (cand_map.count(node_id) == 0) {
                cand_map[node_id] = Candidate(node_id);
            }
            cand_map.at(node_id).rank_n++;
        }

        Vector<int> result(cand_map.size());
        Vector<Candidate> candidates(cand_map.size());
        int idx = 0;
        for (auto & p : cand_map) {
            result[idx] = p.second.id;
            candidates[idx++] = std::move(p.second);
        }

        int sum_dist = 0;
        for (int j = 1; j < candidates.size(); ++j) {
            const auto & b = candidates[j];
            for (int i = 0; i < j; ++i) {
                const auto & a = candidates[i];
                auto pair_dist = dist(a.coord, b.coord);
                // multiply distance by number of rank pairs across the two nodes
                sum_dist += pair_dist * (a.rank_n * b.rank_n);
                if (flag_very_verbose) {
                    Print() << "    Distance from " << a.id
                            << " to " << b.id
                            << ": " << pair_dist << std::endl;
                }
            }
        }
        double score = (nbh_rank_n > 1) ? (static_cast<double>(sum_dist) / pair_n(nbh_rank_n)) : 0;
        return std::make_pair(std::move(result), score);
    }

    // do a local search starting at current node
    std::pair<Vector<int>, double>
    search_local_nbh(int rank_me, const Vector<int> & sg_node_ids, int nbh_rank_n)
    {
        BL_PROFILE("Machine::search_local_nbh()");

        Print() << "Machine::search_local_nbh() called ..." << std::endl;

        Vector<int> result;

        // construct map of node candidates to select
        std::map<int, Candidate> candidates;
        for (auto node_id : sg_node_ids) {
            if (candidates.count(node_id) == 0) {
                candidates[node_id] = Candidate(node_id);
            }
            candidates.at(node_id).rank_n++;
        }

        if (flag_very_verbose) {
            Print() << "  Candidates:" << std::endl;
            for (const auto & p : candidates) {
                const auto & cand = p.second;
                Print() << "    " << cand.id << " : " << to_str(cand.coord)
                        << ": " << cand.rank_n << " ranks" << std::endl;
            }
        }

        AMREX_ASSERT(rank_me >= 0 && rank_me < sg_node_ids.size());
        Candidate cur_node = std::move(candidates.at(sg_node_ids[rank_me]));
        candidates.erase(cur_node.id);

        // add source_node
        result.push_back(cur_node.id);
        int total_rank_n = cur_node.rank_n;
        int total_pairs_dist = 0;
        if (flag_verbose) {
            Print() << "  Added " << cur_node.id
                    << ": " << to_str(cur_node.coord)
                    << ", ranks: " << cur_node.rank_n
                    << ", total ranks: " << total_rank_n
                    << ", avg dist: " << 0 << std::endl;
        }
        if (total_rank_n >= nbh_rank_n) {
            return {std::move(result), 0};
        }

        double min_avg_dist;
        while (total_rank_n < nbh_rank_n)
        {
            min_avg_dist = std::numeric_limits<double>::max();
            Candidate * next_node = nullptr;
            // update candidates with their pairwise rank distances to cur_node
            for (auto & p : candidates) {
                Candidate & cand_node = p.second;
                auto cand_dist = dist(cand_node.coord, cur_node.coord);
                // multiply distance by number of rank pairs across the two nodes
                cand_node.sum_dist += cand_dist * (cand_node.rank_n * cur_node.rank_n);
                double avg_dist = static_cast<double>(cand_node.sum_dist + total_pairs_dist) /
                                  pair_n(cand_node.rank_n + total_rank_n);
                if (flag_very_verbose) {
                    Print() << "    Distance from " << cand_node.id
                            << " to " << cur_node.id
                            << ": " << cand_dist
                            << ", candidate avg: " << avg_dist << std::endl;
                }
                // keep track of what should be the next node to add
                if (avg_dist < min_avg_dist) {
                    next_node = &cand_node;
                    min_avg_dist = avg_dist;
                }
            }
            cur_node = std::move(*next_node);
            next_node = nullptr;
            candidates.erase(cur_node.id);

            // add cur_node to result
            result.push_back(cur_node.id);
            total_rank_n += cur_node.rank_n;
            total_pairs_dist += cur_node.sum_dist;

            if (flag_verbose) {
                Print() << "  Added " << cur_node.id
                        << ": " << to_str(cur_node.coord)
                        << ", ranks: " << cur_node.rank_n
                        << ", total ranks: " << total_rank_n
                        << ", avg dist: " << min_avg_dist << std::endl;
            }
        }

        return std::make_pair(std::move(result), min_avg_dist);
    }
};

std::unique_ptr<Machine> the_machine;

}

namespace amrex {
namespace machine {

void Initialize () {
    the_machine.reset(new Machine());
    amrex::ExecOnFinalize(machine::Finalize);
}

void Finalize () {
    the_machine.reset();
}

Vector<int> find_best_nbh (int rank_n, bool flag_local_ranks) {
    AMREX_ASSERT(the_machine);
    return the_machine->find_best_nbh(rank_n, flag_local_ranks);
}

}}
