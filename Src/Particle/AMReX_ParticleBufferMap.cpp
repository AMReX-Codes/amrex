#include <AMReX_ParticleBufferMap.H>

using namespace amrex;

ParticleBufferMap::ParticleBufferMap (const ParGDBBase* a_gdb)
{
    define(a_gdb);
}

void ParticleBufferMap::define (const ParGDBBase* a_gdb)
{
    BL_PROFILE("ParticleBufferMap::define");

    m_defined = true;

    int num_levels = a_gdb->finestLevel()+1;
    m_ba.resize(0);
    m_dm.resize(0);
    m_ba.resize(num_levels);
    m_dm.resize(num_levels);
    for (int lev = 0; lev < num_levels; ++lev)
    {
        m_ba[lev] = a_gdb->ParticleBoxArray(lev);
        m_dm[lev] = a_gdb->ParticleDistributionMap(lev);
    }

    m_lev_offsets.resize(0);
    m_lev_offsets.push_back(0);
    for (int lev = 0; lev < num_levels; ++lev)
        m_lev_offsets.push_back(m_lev_offsets.back() + m_ba[lev].size());

    int num_buckets = m_lev_offsets.back();

    m_bucket_to_gid.resize(0);
    m_bucket_to_gid.resize(num_buckets);
    m_bucket_to_lev.resize(0);
    m_bucket_to_lev.resize(num_buckets);
    m_bucket_to_pid.resize(0);
    m_bucket_to_pid.resize(num_buckets);

    m_lev_gid_to_bucket.resize(0);
    m_lev_gid_to_bucket.resize(num_buckets);

    using ThreeIntTuple = std::tuple<int, int, int>;
    std::vector<ThreeIntTuple> box_lev_proc_ids;

    for (int lev = 0; lev < num_levels; ++lev) {
        for (int i = 0; i < m_ba[lev].size(); ++i) {
            int rank = ParallelContext::global_to_local_rank(m_dm[lev][i]);
            box_lev_proc_ids.push_back(std::make_tuple(i, lev, rank));
        }
    }

    std::sort(box_lev_proc_ids.begin(), box_lev_proc_ids.end(),
              [](const ThreeIntTuple& a, const ThreeIntTuple& b) -> bool
              {
                  int pid_a = std::get<2>(a);
                  int pid_b = std::get<2>(b);
                  if (pid_a != pid_b) return pid_a < pid_b;

                  int lev_a = std::get<1>(a);
                  int lev_b = std::get<1>(b);
                  if (lev_a != lev_b) return lev_a < lev_b;

                  int gid_a = std::get<0>(a);
                  int gid_b = std::get<0>(b);
                  if (gid_a != gid_b) return gid_a < gid_b;

                  return false;
              });

    int bucket_index = 0;
    for (int lev = 0; lev < num_levels; ++lev) {
        for (int i = 0; i < m_ba[lev].size(); ++i) {
            m_bucket_to_gid[bucket_index] = std::get<0>(box_lev_proc_ids[bucket_index]);
            m_bucket_to_lev[bucket_index] = std::get<1>(box_lev_proc_ids[bucket_index]);
            m_bucket_to_pid[bucket_index] = std::get<2>(box_lev_proc_ids[bucket_index]);
            ++bucket_index;
        }
    }

    m_proc_box_counts.resize(0);
    m_proc_box_counts.resize(ParallelContext::NProcsSub(), 0);

    for (int i = 0; i < num_buckets; ++i)
    {
        int lev = m_bucket_to_lev[i];
        int pid = m_bucket_to_pid[i];
        m_lev_gid_to_bucket[m_lev_offsets[lev] + m_bucket_to_gid[i]] = i;
        m_proc_box_counts[pid]++;
    }

    m_proc_box_offsets.resize(0);
    m_proc_box_offsets.push_back(0);
    for (auto count : m_proc_box_counts)
        m_proc_box_offsets.push_back(m_proc_box_offsets.back() + count);

    d_bucket_to_pid.resize(0);
    d_bucket_to_pid.resize(num_buckets);

    d_lev_gid_to_bucket.resize(0);
    d_lev_gid_to_bucket.resize(num_buckets);

    d_lev_offsets.resize(0);
    d_lev_offsets.resize(m_lev_offsets.size());

    Gpu::copy(Gpu::hostToDevice, m_lev_gid_to_bucket.begin(),m_lev_gid_to_bucket.end(),d_lev_gid_to_bucket.begin());
    Gpu::copy(Gpu::hostToDevice, m_lev_offsets.begin(),m_lev_offsets.end(),d_lev_offsets.begin());
    Gpu::copy(Gpu::hostToDevice, m_bucket_to_pid.begin(),m_bucket_to_pid.end(),d_bucket_to_pid.begin());
}

bool ParticleBufferMap::isValid (const ParGDBBase* a_gdb) const
{
    if (!m_defined) return false;

    int num_levs = a_gdb->finestLevel() + 1;
    if (num_levs != m_ba.size()) return false;

    bool valid = true;
    for (int lev = 0; lev < num_levs; ++lev)
    {
        bool same_ba = BoxArray::SameRefs(a_gdb->ParticleBoxArray(lev), m_ba[lev]);
        bool same_dm = DistributionMapping::SameRefs(a_gdb->ParticleDistributionMap(lev), m_dm[lev]);
        bool level_valid = same_ba && same_dm;
        valid = valid && level_valid;
    }

    return valid;
}
