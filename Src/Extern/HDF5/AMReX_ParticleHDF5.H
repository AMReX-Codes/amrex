#ifndef AMREX_PARTICLEHDF5_H
#define AMREX_PARTICLEHDF5_H
#include <AMReX.H>
#include <AMReX_Config.H>

#include <AMReX_WriteBinaryParticleData.H>

#ifdef AMREX_USE_HDF5
#include "hdf5.h"

#ifdef AMREX_USE_HDF5_ZFP
#include "H5Zzfp_lib.h"
#include "H5Zzfp_props.h"
#endif

#ifdef AMREX_USE_HDF5_SZ
#include "H5Z_SZ.h"
#endif

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CheckpointHDF5 (const std::string& dir,
                  const std::string& name, bool /*is_checkpoint*/,
                  const Vector<std::string>& real_comp_names,
                  const Vector<std::string>& int_comp_names,
                  const std::string& compression) const
{
    Vector<int> write_real_comp;
    Vector<std::string> tmp_real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        if (real_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "real_comp" << i;
            tmp_real_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_real_comp_names.push_back(real_comp_names[i]);
        }
    }

    Vector<int> write_int_comp;
    Vector<std::string> tmp_int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        if (int_comp_names.size() == 0)
        {
            std::stringstream ss;
            ss << "int_comp" << i;
            tmp_int_comp_names.push_back(ss.str());
        }
        else
        {
            tmp_int_comp_names.push_back(int_comp_names[i]);
        }
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          tmp_real_comp_names, tmp_int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p) -> int
                          {
                              return p.id() > 0;
                          }, true);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CheckpointHDF5 (const std::string& dir, const std::string& name,
                  const std::string& compression) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p) -> int
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const std::string& compression) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const Vector<std::string>& real_comp_names,
                     const Vector<std::string>& int_comp_names,
                     const std::string& compression) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NumIntComps() );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const Vector<std::string>& real_comp_names,
                     const std::string& compression) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p)
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir,
                     const std::string& name,
                     const Vector<int>& write_real_comp,
                     const Vector<int>& write_int_comp,
                     const std::string& compression) const
{
    AMREX_ASSERT(write_real_comp.size() == NStructReal + NumRealComps());
    AMREX_ASSERT(write_int_comp.size()  == NStructInt  + NArrayInt );

    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::
WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                   const Vector<int>& write_real_comp,
                   const Vector<int>& write_int_comp,
                   const Vector<std::string>& real_comp_names,
                   const Vector<std::string>&  int_comp_names,
                   const std::string& compression) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression,
                          [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                          {
                              return p.id() > 0;
                          });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>>::value>::type*>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const std::string& compression, F&& f) const
{
    Vector<int> write_real_comp;
    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<int> write_int_comp;
    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        write_int_comp.push_back(1);
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names, compression,
                          std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const Vector<std::string>& real_comp_names,
                     const Vector<std::string>& int_comp_names,
                     const std::string& compression, F&& f) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NArrayInt );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression, std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>>::value>::type*>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                     const Vector<std::string>& real_comp_names,
                     const std::string& compression, F&& f) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression, std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFileHDF5 (const std::string& dir,
                     const std::string& name,
                     const Vector<int>& write_real_comp,
                     const Vector<int>& write_int_comp,
                     const std::string& compression, F&& f) const
{
    AMREX_ASSERT(write_real_comp.size() == NStructReal + NumRealComps());
    AMREX_ASSERT(write_int_comp.size()  == NStructInt  + NumIntComps() );

    Vector<std::string> real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        std::stringstream ss;
        ss << "real_comp" << i;
        real_comp_names.push_back(ss.str());
    }

    Vector<std::string> int_comp_names;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    {
        std::stringstream ss;
        ss << "int_comp" << i;
        int_comp_names.push_back(ss.str());
    }

    WriteHDF5ParticleData(dir, name, write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression, std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::
WritePlotFileHDF5 (const std::string& dir, const std::string& name,
                   const Vector<int>& write_real_comp,
                   const Vector<int>& write_int_comp,
                   const Vector<std::string>& real_comp_names,
                   const Vector<std::string>&  int_comp_names,
                   const std::string& compression, F&& f) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");

    WriteHDF5ParticleData(dir, name,
                          write_real_comp, write_int_comp,
                          real_comp_names, int_comp_names,
                          compression, std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteHDF5ParticleData (const std::string& dir, const std::string& name,
                         const Vector<int>& write_real_comp,
                         const Vector<int>& write_int_comp,
                         const Vector<std::string>& real_comp_names,
                         const Vector<std::string>& int_comp_names,
                         const std::string& compression,
                         F&& f, bool is_checkpoint) const
{
    /* HDF5 async is implemented in WriteHDF5ParticleDataSync, enabled by compile flag */
    /* if (AsyncOut::UseAsyncOut()) { */
    /*     WriteHDF5ParticleDataAsync(*this, dir, name, */
    /*                                write_real_comp, write_int_comp, */
    /*                                real_comp_names, int_comp_names); */
    /* } else */
    /* { */
        WriteHDF5ParticleDataSync(*this, dir, name,
                                  write_real_comp, write_int_comp,
                                  real_comp_names, int_comp_names,
                                  compression,
                                  std::forward<F>(f), is_checkpoint);
    /* } */
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CheckpointPreHDF5 ()
{
    if( ! usePrePost) {
        return;
    }

    BL_PROFILE("ParticleContainer::CheckpointPre()");

    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    Long nparticles = 0;
    Long  maxnextid  = ParticleType::NextID();

    for (int lev = 0; lev < m_particles.size();  lev++) {
        const auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            for (int k = 0; k < aos.numParticles(); ++k) {
                const ParticleType& p = aos[k];
                if (p.id() > 0) {
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
                }
            }
        }
    }
    ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);

    ParticleType::NextID(maxnextid);
    ParallelDescriptor::ReduceLongMax(maxnextid, IOProcNumber);

    nparticlesPrePost = nparticles;
    maxnextidPrePost  = maxnextid;

    nParticlesAtLevelPrePost.clear();
    nParticlesAtLevelPrePost.resize(finestLevel() + 1, 0);
    for(int lev(0); lev <= finestLevel(); ++lev) {
        nParticlesAtLevelPrePost[lev] = NumberOfParticlesAtLevel(lev);
    }

    whichPrePost.clear();
    whichPrePost.resize(finestLevel() + 1);
    countPrePost.clear();
    countPrePost.resize(finestLevel() + 1);
    wherePrePost.clear();
    wherePrePost.resize(finestLevel() + 1);

    filePrefixPrePost.clear();
    filePrefixPrePost.resize(finestLevel() + 1);
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CheckpointPostHDF5 ()
{
    if( ! usePrePost) {
        return;
    }

    BL_PROFILE("ParticleContainer::CheckpointPost()");

    const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();
    std::ofstream HdrFile;
    HdrFile.open(HdrFileNamePrePost.c_str(), std::ios::out | std::ios::app);

    for(int lev(0); lev <= finestLevel(); ++lev) {
        ParallelDescriptor::ReduceIntSum (whichPrePost[lev].dataPtr(), whichPrePost[lev].size(), IOProcNumber);
        ParallelDescriptor::ReduceIntSum (countPrePost[lev].dataPtr(), countPrePost[lev].size(), IOProcNumber);
        ParallelDescriptor::ReduceLongSum(wherePrePost[lev].dataPtr(), wherePrePost[lev].size(), IOProcNumber);


        if(ParallelDescriptor::IOProcessor()) {
            for(int j(0); j < whichPrePost[lev].size(); ++j) {
                HdrFile << whichPrePost[lev][j] << ' ' << countPrePost[lev][j] << ' ' << wherePrePost[lev][j] << '\n';
            }

            const bool gotsome = (nParticlesAtLevelPrePost[lev] > 0);
            if(gotsome && doUnlink) {
//            BL_PROFILE_VAR("PC<NNNN>::Checkpoint:unlink", unlink_post);
                // Unlink any zero-length data files.
                Vector<Long> cnt(nOutFilesPrePost,0);

                for(int i(0), N = countPrePost[lev].size(); i < N; ++i) {
                    cnt[whichPrePost[lev][i]] += countPrePost[lev][i];
                }

                for(int i(0), N = cnt.size(); i < N; ++i) {
                    if(cnt[i] == 0) {
                        std::string FullFileName = NFilesIter::FileName(i, filePrefixPrePost[lev]);
                        FileSystem::Remove(FullFileName);
                    }
                }
            }
        }
    }

    if(ParallelDescriptor::IOProcessor()) {
        HdrFile.flush();
        HdrFile.close();
        if( ! HdrFile.good()) {
            amrex::Abort("ParticleContainer::CheckpointPost(): problem writing HdrFile");
        }
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFilePreHDF5 ()
{
    CheckpointPreHDF5();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFilePostHDF5 ()
{
    CheckpointPostHDF5();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteParticlesHDF5 (int lev, hid_t grp,
                      Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                      const Vector<int>& write_real_comp,
                      const Vector<int>& write_int_comp,
                      const std::string& compression,
                      const Vector<std::map<std::pair<int, int>, IntVector>>& particle_io_flags,
                      bool is_checkpoint) const
{
    BL_PROFILE("ParticleContainer::WriteParticlesHDF5()");

    // For a each grid, the tiles it contains
    std::map<int, Vector<int> > tile_map;

    int ret;
    hid_t dxpl_col, dxpl_ind, dcpl_int, dcpl_real;
    dxpl_col = H5Pcreate(H5P_DATASET_XFER);
    dxpl_ind = H5Pcreate(H5P_DATASET_XFER);
#ifdef BL_USE_MPI
    H5Pset_dxpl_mpio(dxpl_col, H5FD_MPIO_COLLECTIVE);
#endif
    dcpl_int  = H5Pcreate(H5P_DATASET_CREATE);
    dcpl_real = H5Pcreate(H5P_DATASET_CREATE);

#if (defined AMREX_USE_HDF5_ZFP) || (defined AMREX_USE_HDF5_SZ)
    const char *chunk_env = NULL;
    std::string mode_env, value_env;
    double comp_value = -1.0;
    hsize_t chunk_dim = 1024;

    chunk_env = getenv("HDF5_CHUNK_SIZE");
    if (chunk_env != NULL)
        chunk_dim = atoi(chunk_env);

    H5Pset_chunk(dcpl_int, 1, &chunk_dim);
    H5Pset_chunk(dcpl_real, 1, &chunk_dim);

    H5Pset_alloc_time(dcpl_int, H5D_ALLOC_TIME_LATE);
    H5Pset_alloc_time(dcpl_real, H5D_ALLOC_TIME_LATE);

    std::string::size_type pos = compression.find('@');
    if (pos != std::string::npos) {
        mode_env = compression.substr(0, pos);
        value_env = compression.substr(pos+1);
        if (!value_env.empty()) {
            comp_value = atof(value_env.c_str());
        }
    }

#ifdef AMREX_USE_HDF5_ZFP
    pos = compression.find("ZFP");
    if (pos != std::string::npos) {
        ret = H5Z_zfp_initialize();
        if (ret < 0) amrex::Abort("ZFP initialize failed!");
    }
#endif

#ifdef AMREX_USE_HDF5_SZ
    pos = compression.find("SZ");
    if (pos != std::string::npos) {
        ret = H5Z_SZ_Init((char*)value_env.c_str());
        if (ret < 0) {
            std::cout << "SZ config file:" << value_env.c_str() << std::endl;
            amrex::Abort("SZ initialize failed, check SZ config file!");
        }
    }
#endif

    if (!mode_env.empty() && mode_env != "None") {
        if (mode_env == "ZLIB") {
            H5Pset_deflate(dcpl_int, (int)comp_value);
            H5Pset_deflate(dcpl_real, (int)comp_value);
        }
#ifdef AMREX_USE_HDF5_ZFP
        else if (mode_env == "ZFP_RATE") {
            H5Pset_zfp_rate(dcpl_int, comp_value);
            H5Pset_zfp_rate(dcpl_real, comp_value);
        }
        else if (mode_env == "ZFP_PRECISION") {
            H5Pset_zfp_precision(dcpl_int, (unsigned int)comp_value);
            H5Pset_zfp_precision(dcpl_real, (unsigned int)comp_value);
        }
        else if (mode_env == "ZFP_ACCURACY") {
            H5Pset_zfp_accuracy(dcpl_int, comp_value);
            H5Pset_zfp_accuracy(dcpl_real, comp_value);
        }
        else if (mode_env == "ZFP_REVERSIBLE") {
            H5Pset_zfp_reversible(dcpl_int);
            H5Pset_zfp_reversible(dcpl_real);
        }
        else if (mode_env == "ZLIB") {
            H5Pset_deflate(dcpl_int, (int)comp_value);
            H5Pset_deflate(dcpl_real, (int)comp_value);
        }
#endif
        if (ParallelDescriptor::MyProc() == 0) {
            std::cout << "\nHDF5 particle checkpoint using " << mode_env << ", " <<
                value_env << ", " << chunk_dim << std::endl;
        }
    }
#endif


    for (const auto& kv : m_particles[lev])
    {
        const int grid = kv.first.first;
        const int tile = kv.first.second;
        tile_map[grid].push_back(tile);
        const auto& pflags = particle_io_flags[lev].at(kv.first);

        // Only write out valid particles.
        count[grid] += particle_detail::countFlags(pflags);
    }

    MFInfo info;
    info.SetAlloc(false);
    MultiFab state(ParticleBoxArray(lev), ParticleDistributionMap(lev), 1,0,info);

    int            my_mfi_cnt = 0;
    ULong          my_mfi_int_total_size = 0, my_mfi_real_total_size = 0, int_size, real_size;
    Vector<int>    all_mfi_cnt(ParallelDescriptor::NProcs());
    Vector<int>    my_mfi_real_size;
    Vector<int>    my_mfi_int_size;
    Vector<int>    my_nparticles;
    Vector<ULong>  all_mfi_real_total_size(ParallelDescriptor::NProcs());
    Vector<ULong>  all_mfi_int_total_size(ParallelDescriptor::NProcs());
    hid_t          offset_id, offset_space, real_mem_space, real_dset_space, real_dset_id;
    hid_t          int_mem_space, int_dset_id, int_dset_space;
    hsize_t        total_mfi = 0, total_real_size = 0, total_int_size = 0, real_file_offset = 0, int_file_offset = 0;
    hsize_t        my_int_offset, my_int_count, my_real_offset, my_real_count;

    // Get the size for each mf so we know the amount of data from each rank
    for (MFIter mfi(state); mfi.isValid(); ++mfi) {
        const int grid = mfi.index();
        if (count[grid] == 0)
            continue;

        int_size  = count[grid] * (2 + NStructInt + NArrayInt);
        my_mfi_int_size.push_back(int_size);
        my_nparticles.push_back(count[grid]);
        my_mfi_int_total_size += int_size;

        real_size = count[grid] * (AMREX_SPACEDIM + NStructReal + NArrayReal);
        my_mfi_real_size.push_back(real_size);
        my_mfi_real_total_size += real_size;
        my_mfi_cnt++;
    }

    #ifdef BL_USE_MPI
    // Collect the number of mf and total size of mf from each rank
    MPI_Allgather(&my_mfi_cnt, 1, ParallelDescriptor::Mpi_typemap<int>::type(), &(all_mfi_cnt[0]), 1,
                  ParallelDescriptor::Mpi_typemap<int>::type(), ParallelDescriptor::Communicator());

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        total_mfi += all_mfi_cnt[i];

    // Create the int data
    MPI_Allgather(&my_mfi_int_total_size, 1, ParallelDescriptor::Mpi_typemap<ULong>::type(),
                  &(all_mfi_int_total_size[0]), 1, ParallelDescriptor::Mpi_typemap<ULong>::type(), ParallelDescriptor::Communicator());
    #else
    all_mfi_cnt[0] = my_mfi_cnt;
    all_mfi_int_total_size[0] = my_mfi_int_total_size;
    #endif

    int_file_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++)
        int_file_offset += all_mfi_int_total_size[i];
    my_int_offset = int_file_offset;
    my_int_count  = 0;

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        total_int_size += all_mfi_int_total_size[i];

    // SZ int compression seems to have issues at the moment
/* #ifdef AMREX_USE_HDF5_SZ */
/*     if (mode_env == "SZ") { */
/*         size_t cd_nelmts; */
/*         unsigned int* cd_values = NULL; */
/*         unsigned filter_config; */
/*         SZ_metaDataToCdArray(&cd_nelmts, &cd_values, SZ_INT32, 0, 0, 0, 0, total_int_size); */
/*         H5Pset_filter(dcpl_int, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values); */
/*     } */
/* #endif */

    int_dset_space = H5Screate_simple(1, &total_int_size, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
    int_dset_id  = H5Dcreate_async(grp, "data:datatype=0", H5T_NATIVE_INT, int_dset_space, H5P_DEFAULT, dcpl_int, H5P_DEFAULT, es_par_g);
#else
    int_dset_id  = H5Dcreate(grp, "data:datatype=0", H5T_NATIVE_INT, int_dset_space, H5P_DEFAULT, dcpl_int, H5P_DEFAULT);
#endif

    H5Sclose(int_dset_space);

    // Create the real data
    #ifdef BL_USE_MPI
    MPI_Allgather(&my_mfi_real_total_size, 1, ParallelDescriptor::Mpi_typemap<ULong>::type(),
                  &(all_mfi_real_total_size[0]), 1, ParallelDescriptor::Mpi_typemap<ULong>::type(), ParallelDescriptor::Communicator());
    #else
    all_mfi_real_total_size[0] = my_mfi_real_total_size;
    #endif

    for (int i = 0; i < ParallelDescriptor::NProcs(); i++)
        total_real_size += all_mfi_real_total_size[i];

#ifdef AMREX_USE_HDF5_SZ
    if (mode_env == "SZ") {
        size_t cd_nelmts;
        unsigned int* cd_values = NULL;
        unsigned filter_config;
        if (sizeof(typename ParticleType::RealType) == 4) {
            SZ_metaDataToCdArray(&cd_nelmts, &cd_values, SZ_FLOAT, 0, 0, 0, 0, total_real_size);
            H5Pset_filter(dcpl_real, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
        else {
            SZ_metaDataToCdArray(&cd_nelmts, &cd_values, SZ_DOUBLE, 0, 0, 0, 0, total_real_size);
            H5Pset_filter(dcpl_real, H5Z_FILTER_SZ, H5Z_FLAG_MANDATORY, cd_nelmts, cd_values);
        }
    }
#endif

    real_dset_space = H5Screate_simple(1, &total_real_size, NULL);
    if (sizeof(typename ParticleType::RealType) == 4) {
#ifdef AMREX_USE_HDF5_ASYNC
        real_dset_id  = H5Dcreate_async(grp, "data:datatype=1", H5T_NATIVE_FLOAT, real_dset_space,
                                  H5P_DEFAULT, dcpl_real, H5P_DEFAULT, es_par_g);
#else
        real_dset_id  = H5Dcreate(grp, "data:datatype=1", H5T_NATIVE_FLOAT, real_dset_space,
                                  H5P_DEFAULT, dcpl_real, H5P_DEFAULT);
#endif
    }
    else {
#ifdef AMREX_USE_HDF5_ASYNC
        real_dset_id  = H5Dcreate_async(grp, "data:datatype=1", H5T_NATIVE_DOUBLE, real_dset_space,
                                  H5P_DEFAULT, dcpl_real, H5P_DEFAULT, es_par_g);
#else
        real_dset_id  = H5Dcreate(grp, "data:datatype=1", H5T_NATIVE_DOUBLE, real_dset_space,
                                  H5P_DEFAULT, dcpl_real, H5P_DEFAULT);
#endif
    }
    H5Sclose(real_dset_space);

    real_file_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++)
        real_file_offset += all_mfi_real_total_size[i];
    my_real_offset = real_file_offset;
    my_real_count  = 0;

    int max_mfi_count = 0, write_count = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++)
        if (max_mfi_count < all_mfi_cnt[i])
            max_mfi_count = all_mfi_cnt[i];


    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();

        if (count[grid] == 0) continue;

        Vector<int> istuff;
        Vector<ParticleReal> rstuff;
        particle_detail::packIOData(istuff, rstuff, *this, lev, grid,
                                    write_real_comp, write_int_comp,
                                    particle_io_flags, tile_map[grid], count[grid],
                                    is_checkpoint);

        my_int_count = istuff.size();
        int_mem_space = H5Screate_simple(1, &my_int_count, NULL);
        /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_int_offset = " << */
        /*         my_int_offset << ", my_int_count = " << my_int_count << ", total_int_size = " << total_int_size << std::endl; */
        int_dset_space = H5Screate_simple(1, &total_int_size, NULL);
        H5Sselect_hyperslab (int_dset_space, H5S_SELECT_SET, &my_int_offset, NULL, &my_int_count, NULL);

#ifdef AMREX_USE_HDF5_ASYNC
        ret = H5Dwrite_async(int_dset_id, H5T_NATIVE_INT, int_mem_space, int_dset_space, dxpl_col, istuff.dataPtr(), es_par_g);
#else
        ret = H5Dwrite(int_dset_id, H5T_NATIVE_INT, int_mem_space, int_dset_space, dxpl_col, istuff.dataPtr());
#endif
        if (ret < 0) amrex::Abort("H5Dwrite int_dset failed!");

        H5Sclose(int_dset_space);
        H5Sclose(int_mem_space);

        my_int_offset += my_int_count;

        my_real_count = rstuff.size();
        real_mem_space = H5Screate_simple(1, &my_real_count, NULL);
        /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_real_offset = " << */
        /*         my_real_offset << ", my_real_count = " << my_real_count << ", total_real_size = " << total_real_size << '\n'; */
        real_dset_space = H5Screate_simple(1, &total_real_size, NULL);
        H5Sselect_hyperslab (real_dset_space, H5S_SELECT_SET, &my_real_offset, NULL, &my_real_count, NULL);
        if (sizeof(typename ParticleType::RealType) == 4)
#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(real_dset_id, H5T_NATIVE_FLOAT, real_mem_space, real_dset_space, dxpl_col, rstuff.dataPtr(), es_par_g);
#else
            ret = H5Dwrite(real_dset_id, H5T_NATIVE_FLOAT, real_mem_space, real_dset_space, dxpl_col, rstuff.dataPtr());
#endif
        else
#ifdef AMREX_USE_HDF5_ASYNC
            ret = H5Dwrite_async(real_dset_id, H5T_NATIVE_DOUBLE, real_mem_space, real_dset_space, dxpl_col, rstuff.dataPtr(), es_par_g);
#else
            ret = H5Dwrite(real_dset_id, H5T_NATIVE_DOUBLE, real_mem_space, real_dset_space, dxpl_col, rstuff.dataPtr());
#endif

        if (ret < 0) amrex::Abort("H5Dwrite real_dset failed!");

        H5Sclose(real_mem_space);
        H5Sclose(real_dset_space);

        my_real_offset += my_real_count;
        write_count++;

    } // end for (mfi)

    // Dummy writes so that every rank participates to every possible H5Dwrite (collective)
    while (write_count < max_mfi_count) {
        int_dset_space = H5Screate_simple(1, &total_int_size, NULL);
        real_dset_space = H5Screate_simple(1, &total_real_size, NULL);
        H5Sselect_none(int_dset_space);
        H5Sselect_none(real_dset_space);

#ifdef AMREX_USE_HDF5_ASYNC
        H5Dwrite_async(int_dset_id, H5T_NATIVE_INT, int_dset_space, int_dset_space, dxpl_col, NULL, es_par_g);
        if (sizeof(typename ParticleType::RealType) == 4)
            H5Dwrite_async(real_dset_id, H5T_NATIVE_FLOAT, real_dset_space, real_dset_space, dxpl_col, NULL, es_par_g);
        else
            H5Dwrite_async(real_dset_id, H5T_NATIVE_DOUBLE, real_dset_space, real_dset_space, dxpl_col, NULL, es_par_g);
#else
        H5Dwrite(int_dset_id, H5T_NATIVE_INT, int_dset_space, int_dset_space, dxpl_col, NULL);
        if (sizeof(typename ParticleType::RealType) == 4)
            H5Dwrite(real_dset_id, H5T_NATIVE_FLOAT, real_dset_space, real_dset_space, dxpl_col, NULL);
        else
            H5Dwrite(real_dset_id, H5T_NATIVE_DOUBLE, real_dset_space, real_dset_space, dxpl_col, NULL);
#endif

        H5Sclose(int_dset_space);
        H5Sclose(real_dset_space);

        write_count++;
    }

#ifdef AMREX_USE_HDF5_ASYNC
    H5Dclose_async(real_dset_id, es_par_g);
    H5Dclose_async(int_dset_id, es_par_g);
#else
    H5Dclose(real_dset_id);
    H5Dclose(int_dset_id);
#endif

    // Create and write the size dataset
    offset_space = H5Screate_simple(1, &total_mfi, NULL);
#ifdef AMREX_USE_HDF5_ASYNC
    offset_id  = H5Dcreate_async(grp, "nparticles_grid", H5T_NATIVE_INT, offset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT, es_par_g);
#else
    offset_id  = H5Dcreate(grp, "nparticles_grid", H5T_NATIVE_INT, offset_space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
#endif

    my_int_offset = 0;
    for (int i = 0; i < ParallelDescriptor::MyProc(); i++)
        my_int_offset += all_mfi_cnt[i];
    my_int_count  = my_mfi_cnt;
    int_mem_space = H5Screate_simple(1, &my_int_count, NULL);
    /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": my_int_offset = " << */
    /*         my_int_offset << ", my_int_count = " << my_int_count << ", total_mfi = " << total_mfi << '\n'; */
    H5Sselect_hyperslab (offset_space, H5S_SELECT_SET, &my_int_offset, NULL, &my_int_count, NULL);

#ifdef AMREX_USE_HDF5_ASYNC
    ret = H5Dwrite_async(offset_id, H5T_NATIVE_INT, int_mem_space, offset_space, dxpl_col, &(my_nparticles[0]), es_par_g);
#else
    ret = H5Dwrite(offset_id, H5T_NATIVE_INT, int_mem_space, offset_space, dxpl_col, &(my_nparticles[0]));
#endif
    if (ret < 0) amrex::Abort("H5Dwrite offset failed!");

    H5Pclose(dcpl_int);
    H5Pclose(dcpl_real);
    H5Pclose(dxpl_col);
    H5Pclose(dxpl_ind);

    H5Sclose(int_mem_space);
    H5Sclose(offset_space);

#ifdef AMREX_USE_HDF5_ASYNC
    H5Dclose_async(offset_id, es_par_g);
#else
    H5Dclose(offset_id);
#endif

    /* std::cout << "Rank " << ParallelDescriptor::MyProc() << ": done WriteParticlesHDF5" << std::endl; */
    return;
} // End WriteParticlesHDF5

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::RestartHDF5 (const std::string& dir, const std::string& file, bool /*is_checkpoint*/)
{
    RestartHDF5(dir, file);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::RestartHDF5 (const std::string& dir, const std::string& file)
{
    BL_PROFILE("ParticleContainer::RestartHDF5()");
    AMREX_ASSERT(!dir.empty());
    AMREX_ASSERT(!file.empty());

    const auto strttime = amrex::second();

    std::string fullname = dir;
    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';

    fullname += file;
    fullname += ".h5";

    hid_t fid, dset, grp, fapl, attr, atype, dspace, int_dset, real_dset;
    int ret;

    fapl = H5Pcreate (H5P_FILE_ACCESS);
#ifdef BL_USE_MPI
    SetHDF5fapl(fapl, ParallelDescriptor::Communicator());
#else
    SetHDF5fapl(fapl);
#endif

    fid = H5Fopen(fullname.c_str(), H5F_ACC_RDONLY, fapl);
    if (fid < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open file: ");
        msg += fullname;
        amrex::Abort(msg.c_str());
    }

    attr = H5Aopen(fid, "version_name", H5P_DEFAULT);
    if (attr < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open version attribute ");
        amrex::Abort(msg.c_str());
    }

    atype = H5Aget_type(attr);
    if (atype < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to get type of attribute ");
        amrex::Abort(msg.c_str());
    }

    size_t attr_len = H5Tget_size(atype);
    if (attr_len == 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to get size of attribute ");
        amrex::Abort(msg.c_str());
    }

    std::string version;
    version.resize(attr_len+1);

    ret = H5Aread(attr, atype, &version[0]);
    H5Tclose(atype);

    H5Aclose(attr);

    AMREX_ASSERT(!version.empty());

    // What do our version strings mean?
    // "Version_One_Dot_Zero" -- hard-wired to write out in double precision.
    // "Version_One_Dot_One" -- can write out either as either single or double precision.
    // Appended to the latter version string are either "_single" or "_double" to
    // indicate how the particles were written.
    // "Version_Two_Dot_Zero" -- this is the AMReX particle file format
    // "Version_Two_Dot_One" -- expanded particle ids to allow for 2**39-1 per proc
    std::string how;
    bool convert_ids = false;
    if (version.find("Version_Two_Dot_One") != std::string::npos) {
        convert_ids = true;
    }
    if (version.find("_single") != std::string::npos) {
        how = "single";
    }
    else if (version.find("_double") != std::string::npos) {
        how = "double";
    }
    else {
        std::string msg("ParticleContainer::Restart(): bad version string: ");
        msg += version;
        amrex::Error(version.c_str());
    }

    std::string gname = "Chombo_global";
    std::string aname = "SpaceDim";
    int dm;
    grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
    if (grp < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
        msg += gname;
        amrex::Abort(msg.c_str());
    }
    ret = ReadHDF5Attr(grp, aname.c_str(), &dm, H5T_NATIVE_INT);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    H5Gclose(grp);

    aname = "num_component_real";
    int nr;
    ret = ReadHDF5Attr(fid, aname.c_str(), &nr, H5T_NATIVE_INT);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    if (nr != NStructReal + NumRealComps())
        amrex::Abort("ParticleContainer::Restart(): nr != NStructReal + NumRealComps()");

    aname = "num_component_int";
    int ni;
    ret = ReadHDF5Attr(fid, aname.c_str(), &ni, H5T_NATIVE_INT);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    if (ni != NStructInt + NumIntComps())
        amrex::Abort("ParticleContainer::Restart(): ni != NStructInt");

    aname = "nparticles";
    Long nparticles;
    ret = ReadHDF5Attr(fid, aname.c_str(), &nparticles, H5T_NATIVE_LLONG);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(nparticles >= 0);

    aname = "maxnextid";
    Long maxnextid;
    ret = ReadHDF5Attr(fid, aname.c_str(), &maxnextid, H5T_NATIVE_LLONG);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(maxnextid > 0);
    ParticleType::NextID(maxnextid);

    aname = "finest_level";
    int finest_level_in_file;
    ret = ReadHDF5Attr(fid, aname.c_str(), &finest_level_in_file, H5T_NATIVE_INT);
    if (ret < 0) {
        std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
        msg += aname;
        amrex::Abort(msg.c_str());
    }
    AMREX_ASSERT(finest_level_in_file >= 0);

    hid_t comp_dtype = H5Tcreate (H5T_COMPOUND, 2 * AMREX_SPACEDIM * sizeof(int));
    if (1 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 1 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (2 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 3 * sizeof(int), H5T_NATIVE_INT);
    }
    else if (3 == AMREX_SPACEDIM) {
        H5Tinsert (comp_dtype, "lo_i", 0 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_j", 1 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "lo_k", 2 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_i", 3 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_j", 4 * sizeof(int), H5T_NATIVE_INT);
        H5Tinsert (comp_dtype, "hi_k", 5 * sizeof(int), H5T_NATIVE_INT);
    }

    // Determine whether this is a dual-grid restart or not.
    Vector<BoxArray> particle_box_arrays(finest_level_in_file + 1);
    bool dual_grid = false;

    for (int lev = 0; lev <= finest_level_in_file; lev++)
    {
        if (lev > finestLevel())
        {
            dual_grid = true;
            break;
        }

        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }

        dset = H5Dopen(grp, "boxes", H5P_DEFAULT);
        dspace = H5Dget_space(dset);
        hsize_t ngrid;
        H5Sget_simple_extent_dims(dspace, &ngrid, NULL);

        Vector<int> boxes(ngrid*AMREX_SPACEDIM*2);
        ret = H5Dread(dset, comp_dtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(boxes[0]));
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read nparticles_grid dataset");
            amrex::Abort(msg.c_str());
        }
        H5Sclose(dspace);
        H5Dclose(dset);
        H5Gclose(grp);

        /* particle_box_arrays[lev].readFrom(phdr_file); */
        particle_box_arrays[lev] = ParticleBoxArray(lev);
        for (int i = 0; i < (int)ngrid; i++) {
            Box tmp ({AMREX_D_DECL(boxes[i*AMREX_SPACEDIM], boxes[i*AMREX_SPACEDIM+1], boxes[i*AMREX_SPACEDIM+2])}, {AMREX_D_DECL(boxes[i*AMREX_SPACEDIM*2], boxes[i*AMREX_SPACEDIM*2+1], boxes[i*AMREX_SPACEDIM*2+2])});
            particle_box_arrays[lev][i] = tmp;
        }

        if (! particle_box_arrays[lev].CellEqual(ParticleBoxArray(lev)))
            dual_grid = true;
    }

    if (dual_grid) {
        for (int lev = 0; lev <= finestLevel(); lev++) {
            SetParticleBoxArray(lev, particle_box_arrays[lev]);
            DistributionMapping pdm(particle_box_arrays[lev]);
            SetParticleDistributionMap(lev, pdm);
        }
    }

    Vector<int> ngrids(finest_level_in_file+1);
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }

        aname = "ngrids";
        ret = ReadHDF5Attr(grp, aname.c_str(), &ngrids[lev], H5T_NATIVE_INT);
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read attribute ");
            msg += aname;
            amrex::Abort(msg.c_str());
        }

        AMREX_ASSERT(ngrids[lev] > 0);
        if (lev <= finestLevel()) {
            AMREX_ASSERT(ngrids[lev] == int(ParticleBoxArray(lev).size()));
        }

        H5Gclose(grp);
    }

    resizeData();

    if (finest_level_in_file > finestLevel()) {
        m_particles.resize(finest_level_in_file+1);
    }

    for (int lev = 0; lev <= finest_level_in_file; lev++) {

        gname = "level_" + std::to_string(lev);
        grp = H5Gopen(fid, gname.c_str(), H5P_DEFAULT);
        if (grp < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open group ");
            msg += gname;
            amrex::Abort(msg.c_str());
        }

        dset = H5Dopen(grp, "nparticles_grid", H5P_DEFAULT);
        Vector<int> count(ngrids[lev]);
        ret = H5Dread(dset, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &(count[0]));
        if (ret < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to read nparticles_grid dataset");
            amrex::Abort(msg.c_str());
        }
        H5Dclose(dset);

        Vector<hsize_t>  offset(ngrids[lev]);
        offset[0] = 0;
        for (int i = 1; i < ngrids[lev]; i++) {
            offset[i] = offset[i-1] + count[i-1];
        }

        int_dset  = H5Dopen(grp, "data:datatype=0", H5P_DEFAULT);
        if (int_dset < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open int dataset");
            amrex::Abort(msg.c_str());
        }
        real_dset = H5Dopen(grp, "data:datatype=1", H5P_DEFAULT);
        if (real_dset < 0) {
            std::string msg("ParticleContainer::RestartHDF5(): unable to open int dataset");
            amrex::Abort(msg.c_str());
        }

        Vector<int> grids_to_read;
        if (lev <= finestLevel()) {
            for (MFIter mfi(*m_dummy_mf[lev]); mfi.isValid(); ++mfi) {
                grids_to_read.push_back(mfi.index());
            }
        } else {

            // we lost a level on restart. we still need to read in particles
            // on finer levels, and put them in the right place via Redistribute()

            const int rank = ParallelDescriptor::MyProc();
            const int NReaders = MaxReaders();
            if (rank >= NReaders) return;

            const int Navg = ngrids[lev] / NReaders;
            const int Nleft = ngrids[lev] - Navg * NReaders;

            int lo, hi;
            if (rank < Nleft) {
                lo = rank*(Navg + 1);
                hi = lo + Navg + 1;
            }
            else {
                lo = rank * Navg + Nleft;
                hi = lo + Navg;
            }

            for (int i = lo; i < hi; ++i) {
                grids_to_read.push_back(i);
            }
        }

        for(int igrid = 0; igrid < static_cast<int>(grids_to_read.size()); ++igrid) {
            const int grid = grids_to_read[igrid];

            if (how == "single") {
                ReadParticlesHDF5<float>(offset[grid], count[grid], grid, lev, int_dset, real_dset, finest_level_in_file, convert_ids);
            }
            else if (how == "double") {
                ReadParticlesHDF5<double>(offset[grid], count[grid], grid, lev, int_dset, real_dset, finest_level_in_file, convert_ids);
            }
            else {
                std::string msg("ParticleContainer::Restart(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }
        }

        H5Dclose(int_dset);
        H5Dclose(real_dset);
        H5Gclose(grp);
    } // end for level

    H5Fclose(fid);

    Redistribute();

    AMREX_ASSERT(OK());

    if (m_verbose > 1) {
        auto stoptime = amrex::second() - strttime;
        ParallelDescriptor::ReduceRealMax(stoptime, ParallelDescriptor::IOProcessorNumber());
        amrex::Print() << "ParticleContainer::Restart() time: " << stoptime << '\n';
    }
}

// Read a batch of particles from the checkpoint file
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class RTYPE>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::ReadParticlesHDF5 (hsize_t offset, hsize_t cnt, int grd, int lev,
                     hid_t int_dset, hid_t real_dset, int finest_level_in_file,
                     bool convert_ids)
{
    BL_PROFILE("ParticleContainer::ReadParticlesHDF5()");
    AMREX_ASSERT(cnt > 0);
    AMREX_ASSERT(lev < int(m_particles.size()));

    hid_t int_dspace, int_fspace, real_dspace, real_fspace;

    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    const int iChunkSize = 2 + NStructInt + NumIntComps();
    Vector<int> istuff(cnt*iChunkSize);
    int_fspace = H5Dget_space(int_dset);
    hsize_t int_cnt = cnt*iChunkSize;
    hsize_t int_offset = offset*iChunkSize;
    int_dspace = H5Screate_simple(1, &int_cnt, NULL);
    H5Sselect_hyperslab (int_fspace, H5S_SELECT_SET, &int_offset, NULL, &int_cnt, NULL);
    H5Dread(int_dset, H5T_NATIVE_INT, int_dspace, int_fspace, H5P_DEFAULT, istuff.dataPtr());

    H5Sclose(int_fspace);
    H5Sclose(int_dspace);

    // Then the real data in binary.
    const int rChunkSize = AMREX_SPACEDIM + NStructReal + NumRealComps();
    Vector<RTYPE> rstuff(cnt*rChunkSize);
    real_fspace = H5Dget_space(real_dset);
    hsize_t real_cnt = cnt*rChunkSize;
    hsize_t real_offset = offset*rChunkSize;
    real_dspace = H5Screate_simple(1, &real_cnt, NULL);
    H5Sselect_hyperslab (real_fspace, H5S_SELECT_SET, &real_offset, NULL, &real_cnt, NULL);
    if (sizeof(RTYPE) == 4)
        H5Dread(real_dset, H5T_NATIVE_FLOAT, real_dspace, real_fspace, H5P_DEFAULT, rstuff.dataPtr());
    else
        H5Dread(real_dset, H5T_NATIVE_DOUBLE, real_dspace, real_fspace, H5P_DEFAULT, rstuff.dataPtr());

    H5Sclose(real_fspace);
    H5Sclose(real_dspace);

    // Now reassemble the particles.
    int*   iptr = istuff.dataPtr();
    RTYPE* rptr = rstuff.dataPtr();

    ParticleType p;
    ParticleLocData pld;

    Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
    host_particles.reserve(15);
    host_particles.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<ParticleReal> > > > host_real_attribs;
    host_real_attribs.reserve(15);
    host_real_attribs.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<int> > > > host_int_attribs;
    host_int_attribs.reserve(15);
    host_int_attribs.resize(finest_level_in_file+1);

    for (hsize_t i = 0; i < cnt; i++) {
        if (convert_ids) {
            std::int32_t  xi, yi;
            std::uint32_t xu, yu;
            xi = iptr[0];
            yi = iptr[1];
            std::memcpy(&xu, &xi, sizeof(xi));
            std::memcpy(&yu, &yi, sizeof(yi));
            p.m_idcpu = ((std::uint64_t)xu) << 32 | yu;
        } else {
            p.id()   = iptr[0];
            p.cpu()  = iptr[1];
        }
        iptr += 2;

        for (int j = 0; j < NStructInt; j++)
        {
            p.idata(j) = *iptr;
            ++iptr;
        }

        AMREX_ASSERT(p.id() > 0);

        AMREX_D_TERM(p.pos(0) = ParticleReal(rptr[0]);,
                     p.pos(1) = ParticleReal(rptr[1]);,
                     p.pos(2) = ParticleReal(rptr[2]););

        rptr += AMREX_SPACEDIM;

        for (int j = 0; j < NStructReal; j++)
        {
            p.rdata(j) = ParticleReal(*rptr);
            ++rptr;
        }

        locateParticle(p, pld, 0, finestLevel(), 0);

        std::pair<int, int> ind(grd, pld.m_tile);

        host_real_attribs[lev][ind].resize(NumRealComps());
        host_int_attribs[lev][ind].resize(NumIntComps());

        // add the struct
        host_particles[lev][ind].push_back(p);

        // add the real...
        for (int icomp = 0; icomp < NumRealComps(); icomp++) {
            host_real_attribs[lev][ind][icomp].push_back(ParticleReal(*rptr));
            ++rptr;
        }

        // ... and int array data
        for (int icomp = 0; icomp < NumIntComps(); icomp++) {
            host_int_attribs[lev][ind][icomp].push_back(*iptr);
            ++iptr;
        }
    }

    for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
      {
        for (auto& kv : host_particles[host_lev]) {
          auto grid = kv.first.first;
          auto tile = kv.first.second;
          const auto& src_tile = kv.second;

          auto& dst_tile = DefineAndReturnParticleTile(host_lev, grid, tile);
          auto old_size = dst_tile.GetArrayOfStructs().size();
          auto new_size = old_size + src_tile.size();
          dst_tile.resize(new_size);

          Gpu::copyAsync(Gpu::hostToDevice, src_tile.begin(), src_tile.end(),
                         dst_tile.GetArrayOfStructs().begin() + old_size);

          for (int i = 0; i < NumRealComps(); ++i) {
              Gpu::copyAsync(Gpu::hostToDevice,
                             host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                             host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                             dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
          }

          for (int i = 0; i < NumIntComps(); ++i) {
              Gpu::copyAsync(Gpu::hostToDevice,
                             host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
                             host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
                             dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
          }
        }
      }

    Gpu::streamSynchronize();
}

#endif /*AMREX_USE_HDF5*/
#endif /*AMREX_PARTICLEHDF5_H*/
