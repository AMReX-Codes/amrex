#include <AMReX_ParticleBufferMap.H>
#include <AMReX_ParticleUtil.H>

namespace amrex {

ParticleBufferMap::ParticleBufferMap (const ParGDBBase* a_gdb)
{
    define(a_gdb, false, IntVect(AMREX_D_DECL(1024000, 1024000, 1024000)));
}

ParticleBufferMap::ParticleBufferMap (const ParGDBBase* a_gdb, bool a_do_tiling,
                                      const IntVect& a_tile_size)
{
    define(a_gdb, a_do_tiling, a_tile_size);
}

void ParticleBufferMap::define (const ParGDBBase* a_gdb)
{
    define(a_gdb, false, IntVect(AMREX_D_DECL(1024000, 1024000, 1024000)));
}

void ParticleBufferMap::define (const ParGDBBase* a_gdb, bool a_do_tiling,
                                const IntVect& a_tile_size)
{
    BL_PROFILE("ParticleBufferMap::define");

    m_defined = true;
    m_do_tiling = a_do_tiling;
    m_tile_size = a_tile_size;

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
    for (int lev = 0; lev < num_levels; ++lev) {
        m_lev_offsets.push_back(m_lev_offsets.back() + static_cast<int>(m_ba[lev].size()) + 1);
    }

    m_gid_offsets.resize(m_lev_offsets.back(), 0);

    int num_buckets = 0;
    for (int lev = 0; lev < num_levels; ++lev) {
        int level_offset = m_lev_offsets[lev];
        for (int gid = 0; gid < m_ba[lev].size(); ++gid) {
            m_gid_offsets[level_offset + gid] = num_buckets;
            num_buckets += numTilesInBox(m_ba[lev][gid], m_do_tiling, m_tile_size);
        }
        m_gid_offsets[level_offset + static_cast<int>(m_ba[lev].size())] = num_buckets;
    }

    m_bucket_to_gid.resize(0);
    m_bucket_to_gid.resize(num_buckets);
    m_bucket_to_tid.resize(0);
    m_bucket_to_tid.resize(num_buckets);
    m_bucket_to_lev.resize(0);
    m_bucket_to_lev.resize(num_buckets);
    m_bucket_to_pid.resize(0);
    m_bucket_to_pid.resize(num_buckets);

    m_lev_gid_tid_to_bucket.resize(0);
    m_lev_gid_tid_to_bucket.resize(num_buckets);

    using FourIntTuple = std::tuple<int, int, int, int>;
    std::vector<FourIntTuple> box_tile_lev_proc_ids;

    for (int lev = 0; lev < num_levels; ++lev) {
        for (int i = 0; i < m_ba[lev].size(); ++i) {
            int rank = ParallelContext::global_to_local_rank(m_dm[lev][i]);
            int num_tiles = numTilesInBox(m_ba[lev][i], m_do_tiling, m_tile_size);
            for (int tid = 0; tid < num_tiles; ++tid) {
                box_tile_lev_proc_ids.emplace_back(i, tid, lev, rank);
            }
        }
    }

    std::sort(box_tile_lev_proc_ids.begin(), box_tile_lev_proc_ids.end(),
              [](const FourIntTuple& a, const FourIntTuple& b) -> bool
              {
                  int pid_a = std::get<3>(a);
                  int pid_b = std::get<3>(b);
                  if (pid_a != pid_b) { return pid_a < pid_b; }

                  int lev_a = std::get<2>(a);
                  int lev_b = std::get<2>(b);
                  if (lev_a != lev_b) { return lev_a < lev_b; }

                  int gid_a = std::get<0>(a);
                  int gid_b = std::get<0>(b);
                  if (gid_a != gid_b) { return gid_a < gid_b; }

                  int tid_a = std::get<1>(a);
                  int tid_b = std::get<1>(b);
                  if (tid_a != tid_b) { return tid_a < tid_b; }

                  return false;
              });

    int bucket_index = 0;
    for (int i = 0; i < num_buckets; ++i) {
        m_bucket_to_gid[bucket_index] = std::get<0>(box_tile_lev_proc_ids[bucket_index]);
        m_bucket_to_tid[bucket_index] = std::get<1>(box_tile_lev_proc_ids[bucket_index]);
        m_bucket_to_lev[bucket_index] = std::get<2>(box_tile_lev_proc_ids[bucket_index]);
        m_bucket_to_pid[bucket_index] = std::get<3>(box_tile_lev_proc_ids[bucket_index]);
        ++bucket_index;
    }

    m_proc_box_counts.resize(0);
    m_proc_box_counts.resize(ParallelContext::NProcsSub(), 0);

    for (int i = 0; i < num_buckets; ++i)
    {
        int lev = m_bucket_to_lev[i];
        int gid = m_bucket_to_gid[i];
        int tid = m_bucket_to_tid[i];
        int pid = m_bucket_to_pid[i];
        m_lev_gid_tid_to_bucket[m_gid_offsets[m_lev_offsets[lev] + gid] + tid] = i;
        m_proc_box_counts[pid]++;
    }

    m_proc_box_offsets.resize(0);
    m_proc_box_offsets.push_back(0);
    for (auto count : m_proc_box_counts) {
        m_proc_box_offsets.push_back(m_proc_box_offsets.back() + count);
    }

    d_bucket_to_pid.resize(0);
    d_bucket_to_pid.resize(num_buckets);

    d_lev_gid_tid_to_bucket.resize(0);
    d_lev_gid_tid_to_bucket.resize(num_buckets);

    d_lev_offsets.resize(0);
    d_lev_offsets.resize(m_lev_offsets.size());

    d_gid_offsets.resize(0);
    d_gid_offsets.resize(m_gid_offsets.size());

    Gpu::copyAsync(Gpu::hostToDevice, m_lev_gid_tid_to_bucket.begin(),m_lev_gid_tid_to_bucket.end(),d_lev_gid_tid_to_bucket.begin());
    Gpu::copyAsync(Gpu::hostToDevice, m_lev_offsets.begin(),m_lev_offsets.end(),d_lev_offsets.begin());
    Gpu::copyAsync(Gpu::hostToDevice, m_gid_offsets.begin(),m_gid_offsets.end(),d_gid_offsets.begin());
    Gpu::copyAsync(Gpu::hostToDevice, m_bucket_to_pid.begin(),m_bucket_to_pid.end(),d_bucket_to_pid.begin());
    Gpu::streamSynchronize();
}

bool ParticleBufferMap::isValid (const ParGDBBase* a_gdb) const
{
    return isValid(a_gdb, false, IntVect(AMREX_D_DECL(1024000, 1024000, 1024000)));
}

bool ParticleBufferMap::isValid (const ParGDBBase* a_gdb, bool a_do_tiling,
                                 const IntVect& a_tile_size) const
{
    if (!m_defined) { return false; }
    if (m_do_tiling != a_do_tiling) { return false; }
    if (m_tile_size != a_tile_size) { return false; }

    int num_levs = a_gdb->finestLevel() + 1;
    if (num_levs != m_ba.size()) { return false; }

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

}
