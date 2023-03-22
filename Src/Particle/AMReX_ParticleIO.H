#ifndef AMREX_PARTICLEIO_H
#define AMREX_PARTICLEIO_H
#include <AMReX_Config.H>

#include <AMReX_WriteBinaryParticleData.H>

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteParticleRealData (void* data, size_t size, std::ostream& os) const
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        writeFloatData((float*) data, size, os, ParticleRealDescriptor);
    }
    else if (sizeof(typename ParticleType::RealType) == 8) {
        writeDoubleData((double*) data, size, os, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::ReadParticleRealData (void* data, size_t size, std::istream& is)
{
    if (sizeof(typename ParticleType::RealType) == 4) {
        readFloatData((float*) data, size, is, ParticleRealDescriptor);
    }
    else if (sizeof(typename ParticleType::RealType) == 8) {
        readDoubleData((double*) data, size, is, ParticleRealDescriptor);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::Checkpoint (const std::string& dir,
              const std::string& name, bool /*is_checkpoint*/,
              const Vector<std::string>& real_comp_names,
              const Vector<std::string>& int_comp_names) const
{
    Vector<int> write_real_comp;
    Vector<std::string> tmp_real_comp_names;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    {
        write_real_comp.push_back(1);
        if (real_comp_names.empty())
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
        if (int_comp_names.empty())
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

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            tmp_real_comp_names, tmp_int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p) -> int
                            {
                                return p.id() > 0;
                            }, true);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name) const
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

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names,
                 const Vector<std::string>& int_comp_names) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NumIntComps() );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names) const
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

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const ParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir,
                 const std::string& name,
                 const Vector<int>& write_real_comp,
                 const Vector<int>& write_int_comp) const
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

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::
WritePlotFile (const std::string& dir, const std::string& name,
               const Vector<int>& write_real_comp,
               const Vector<int>& write_int_comp,
               const Vector<std::string>& real_comp_names,
               const Vector<std::string>&  int_comp_names) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            [=] AMREX_GPU_HOST_DEVICE (const SuperParticleType& p)
                            {
                                return p.id() > 0;
                            });
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>&>::value>::type*>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name, F&& f) const
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

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names,
                 const Vector<std::string>& int_comp_names, F&& f) const
{
    AMREX_ASSERT(real_comp_names.size() == NStructReal + NumRealComps());
    AMREX_ASSERT( int_comp_names.size() == NStructInt  + NArrayInt );

    Vector<int> write_real_comp;
    for (int i = 0; i < NStructReal + NumRealComps(); ++i) write_real_comp.push_back(1);

    Vector<int> write_int_comp;
    for (int i = 0; i < NStructInt + NumIntComps(); ++i) write_int_comp.push_back(1);

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F, typename std::enable_if<!std::is_same<F, Vector<std::string>>::value>::type*>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir, const std::string& name,
                 const Vector<std::string>& real_comp_names, F&& f) const
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

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFile (const std::string& dir,
                 const std::string& name,
                 const Vector<int>& write_real_comp,
                 const Vector<int>& write_int_comp, F&& f) const
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

    WriteBinaryParticleData(dir, name, write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>::
WritePlotFile (const std::string& dir, const std::string& name,
               const Vector<int>& write_real_comp,
               const Vector<int>& write_int_comp,
               const Vector<std::string>& real_comp_names,
               const Vector<std::string>&  int_comp_names,
               F&& f) const
{
    BL_PROFILE("ParticleContainer::WritePlotFile()");

    WriteBinaryParticleData(dir, name,
                            write_real_comp, write_int_comp,
                            real_comp_names, int_comp_names,
                            std::forward<F>(f));
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
template <class F>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteBinaryParticleData (const std::string& dir, const std::string& name,
                           const Vector<int>& write_real_comp,
                           const Vector<int>& write_int_comp,
                           const Vector<std::string>& real_comp_names,
                           const Vector<std::string>& int_comp_names,
                           F&& f, bool is_checkpoint) const
{
    if (AsyncOut::UseAsyncOut()) {
        WriteBinaryParticleDataAsync(*this, dir, name,
                                     write_real_comp, write_int_comp,
                                     real_comp_names, int_comp_names, is_checkpoint);
    } else
    {
        WriteBinaryParticleDataSync(*this, dir, name,
                                    write_real_comp, write_int_comp,
                                    real_comp_names, int_comp_names,
                                    std::forward<F>(f), is_checkpoint);
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::CheckpointPre ()
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
::CheckpointPost ()
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
::WritePlotFilePre ()
{
    CheckpointPre();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WritePlotFilePost ()
{
    CheckpointPost();
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteParticles (int lev, std::ofstream& ofs, int fnum,
                  Vector<int>& which, Vector<int>& count, Vector<Long>& where,
                  const Vector<int>& write_real_comp,
                  const Vector<int>& write_int_comp,
                  const Vector<std::map<std::pair<int, int>, IntVector>>& particle_io_flags,
                  bool is_checkpoint) const
{
    BL_PROFILE("ParticleContainer::WriteParticles()");

    // For a each grid, the tiles it contains
    std::map<int, Vector<int> > tile_map;

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

    for (MFIter mfi(state); mfi.isValid(); ++mfi)
    {
        const int grid = mfi.index();

        which[grid] = fnum;
        where[grid] = VisMF::FileOffset(ofs);

        if (count[grid] == 0) continue;

        Vector<int> istuff;
        Vector<ParticleReal> rstuff;
        particle_detail::packIOData(istuff, rstuff, *this, lev, grid,
                                    write_real_comp, write_int_comp,
                                    particle_io_flags, tile_map[grid], count[grid], is_checkpoint);

        writeIntData(istuff.dataPtr(), istuff.size(), ofs);
        ofs.flush();  // Some systems require this flush() (probably due to a bug)

        WriteParticleRealData(rstuff.dataPtr(), rstuff.size(), ofs);
        ofs.flush();  // Some systems require this flush() (probably due to a bug)
    }
}


template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::Restart (const std::string& dir, const std::string& file, bool /*is_checkpoint*/)
{
    Restart(dir, file);
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::Restart (const std::string& dir, const std::string& file)
{
    BL_PROFILE("ParticleContainer::Restart()");
    AMREX_ASSERT(!dir.empty());
    AMREX_ASSERT(!file.empty());

    const auto strttime = amrex::second();

    int DATA_Digits_Read(5);
    ParmParse pp("particles");
    pp.queryAdd("datadigits_read",DATA_Digits_Read);

    std::string fullname = dir;
    if (!fullname.empty() && fullname[fullname.size()-1] != '/')
        fullname += '/';
    fullname += file;
    std::string HdrFileName = fullname;
    if (!HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
        HdrFileName += '/';
    HdrFileName += "Header";

    Vector<char> fileCharPtr;
    ParallelDescriptor::ReadAndBcastFile(HdrFileName, fileCharPtr);
    std::string fileCharPtrString(fileCharPtr.dataPtr());
    std::istringstream HdrFile(fileCharPtrString, std::istringstream::in);

    std::string version;
    HdrFile >> version;
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
    if (version.find("Version_One_Dot_Zero") != std::string::npos) {
        how = "double";
    }
    else if (version.find("Version_One_Dot_One")  != std::string::npos ||
             version.find("Version_Two_Dot_Zero") != std::string::npos ||
             version.find("Version_Two_Dot_One") != std::string::npos) {
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
    }
    else {
        std::string msg("ParticleContainer::Restart(): unknown version string: ");
        msg += version;
        amrex::Abort(msg.c_str());
    }

    int dm;
    HdrFile >> dm;
    if (dm != AMREX_SPACEDIM)
        amrex::Abort("ParticleContainer::Restart(): dm != AMREX_SPACEDIM");

    int nr;
    HdrFile >> nr;
    if (nr != NStructReal + NumRealComps())
        amrex::Abort("ParticleContainer::Restart(): nr != NStructReal + NumRealComps()");

    std::string comp_name;
    for (int i = 0; i < nr; ++i)
        HdrFile >> comp_name;

    int ni;
    HdrFile >> ni;
    if (ni != NStructInt + NumIntComps())
        amrex::Abort("ParticleContainer::Restart(): ni != NStructInt");

    for (int i = 0; i < ni; ++i)
        HdrFile >> comp_name;

    bool checkpoint;
    HdrFile >> checkpoint;

    Long nparticles;
    HdrFile >> nparticles;
    AMREX_ASSERT(nparticles >= 0);

    Long maxnextid;
    HdrFile >> maxnextid;
    AMREX_ASSERT(maxnextid > 0);
    ParticleType::NextID(maxnextid);

    int finest_level_in_file;
    HdrFile >> finest_level_in_file;
    AMREX_ASSERT(finest_level_in_file >= 0);

    // Determine whether this is a dual-grid restart or not.
    Vector<DistributionMapping> old_dms(finest_level_in_file + 1);
    Vector<BoxArray> old_bas(finest_level_in_file + 1);
    Vector<BoxArray> particle_box_arrays(finest_level_in_file + 1);
    bool dual_grid = false;

    bool have_pheaders = false;
    for (int lev = 0; lev <= finest_level_in_file; lev++)
    {
        std::string phdr_name = fullname;
        phdr_name = amrex::Concatenate(phdr_name + "/Level_", lev, 1);
        phdr_name += "/Particle_H";

        if (amrex::FileExists(phdr_name)) {
            have_pheaders = true;
            break;
        }
    }

    if (have_pheaders)
    {
        for (int lev = 0; lev <= finest_level_in_file; lev++)
        {
            old_dms[lev] = ParticleDistributionMap(lev);
            old_bas[lev] = ParticleBoxArray(lev);
            std::string phdr_name = fullname;
            phdr_name = amrex::Concatenate(phdr_name + "/Level_", lev, 1);
            phdr_name += "/Particle_H";

            if (! amrex::FileExists(phdr_name)) continue;

            Vector<char> phdr_chars;
            ParallelDescriptor::ReadAndBcastFile(phdr_name, phdr_chars);
            std::string phdr_string(phdr_chars.dataPtr());
            std::istringstream phdr_file(phdr_string, std::istringstream::in);

            if (lev > finestLevel())
            {
                dual_grid = true;
                break;
            }

            particle_box_arrays[lev].readFrom(phdr_file);
            if (! particle_box_arrays[lev].CellEqual(ParticleBoxArray(lev))) dual_grid = true;
        }
    } else // if no particle box array information exists in the file, we assume a single grid restart
    {
        dual_grid = false;
    }

    if (dual_grid) {
        for (int lev = 0; lev <= finestLevel(); lev++) {
            // this can happen if there are no particles at a given level in the checkpoint
            if (particle_box_arrays[lev].empty()) {
                particle_box_arrays[lev] = BoxArray(Geom(lev).Domain());
            }
            SetParticleBoxArray(lev, particle_box_arrays[lev]);
            DistributionMapping pdm(particle_box_arrays[lev]);
            SetParticleDistributionMap(lev, pdm);
        }
    }

    Vector<int> ngrids(finest_level_in_file+1);
    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        HdrFile >> ngrids[lev];
        AMREX_ASSERT(ngrids[lev] > 0);
    }

    resizeData();

    if (finest_level_in_file > finestLevel()) {
        m_particles.resize(finest_level_in_file+1);
    }

    for (int lev = 0; lev <= finest_level_in_file; lev++) {
        Vector<int>  which(ngrids[lev]);
        Vector<int>  count(ngrids[lev]);
        Vector<Long> where(ngrids[lev]);
        for (int i = 0; i < ngrids[lev]; i++) {
            HdrFile >> which[i] >> count[i] >> where[i];
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

            if (count[grid] <= 0) continue;

            // The file names in the header file are relative.
            std::string name = fullname;

            if (!name.empty() && name[name.size()-1] != '/')
                name += '/';

            name += "Level_";
            name += amrex::Concatenate("", lev, 1);
            name += '/';
            name += DataPrefix();
            name += amrex::Concatenate("", which[grid], DATA_Digits_Read);

            std::ifstream ParticleFile;

            ParticleFile.open(name.c_str(), std::ios::in | std::ios::binary);

            if (!ParticleFile.good())
                amrex::FileOpenFailed(name);

            ParticleFile.seekg(where[grid], std::ios::beg);

            // Use if constexpr to avoid instantiating the mis-matched
            // type case and triggering the static_assert on the
            // underlying copy calls
            if (how == "single") {
                if constexpr (std::is_same_v<ParticleReal, float>) {
                    ReadParticles<float>(count[grid], grid, lev, ParticleFile, finest_level_in_file, convert_ids);
                } else {
                    amrex::Error("File contains single-precision data, while AMReX is compiled with ParticleReal==double");
                }
            }
            else if (how == "double") {
                if constexpr (std::is_same_v<ParticleReal, double>) {
                    ReadParticles<double>(count[grid], grid, lev, ParticleFile, finest_level_in_file, convert_ids);
                } else {
                    amrex::Error("File contains double-precision data, while AMReX is compiled with ParticleReal==float");
                }
            }
            else {
                std::string msg("ParticleContainer::Restart(): bad parameter: ");
                msg += how;
                amrex::Error(msg.c_str());
            }

            ParticleFile.close();

            if (!ParticleFile.good())
                amrex::Abort("ParticleContainer::Restart(): problem reading particles");
        }
    }

    if (dual_grid) {
        for (int lev = 0; lev <= finest_level_in_file; lev++) {
            SetParticleBoxArray(lev, old_bas[lev]);
            SetParticleDistributionMap(lev, old_dms[lev]);
        }
    }

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
::ReadParticles (int cnt, int grd, int lev, std::ifstream& ifs,
                 int finest_level_in_file, bool convert_ids)
{
    BL_PROFILE("ParticleContainer::ReadParticles()");
    AMREX_ASSERT(cnt > 0);
    AMREX_ASSERT(lev < int(m_particles.size()));

    // First read in the integer data in binary.  We do not store
    // the m_lev and m_grid data on disk.  We can easily recreate
    // that given the structure of the checkpoint file.
    const int iChunkSize = 2 + NStructInt + NumIntComps();
    Vector<int> istuff(cnt*iChunkSize);
    readIntData(istuff.dataPtr(), istuff.size(), ifs, FPC::NativeIntDescriptor());

    // Then the real data in binary.
    const int rChunkSize = AMREX_SPACEDIM + NStructReal + NumRealComps();
    Vector<RTYPE> rstuff(cnt*rChunkSize);
    ReadParticleRealData(rstuff.dataPtr(), rstuff.size(), ifs);

    // Now reassemble the particles.
    int*   iptr = istuff.dataPtr();
    RTYPE* rptr = rstuff.dataPtr();

    ParticleType p;
    ParticleLocData pld;

    Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
    host_particles.reserve(15);
    host_particles.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<RTYPE> > > > host_real_attribs;
    host_real_attribs.reserve(15);
    host_real_attribs.resize(finest_level_in_file+1);

    Vector<std::map<std::pair<int, int>,
                    std::vector<Gpu::HostVector<int> > > > host_int_attribs;
    host_int_attribs.reserve(15);
    host_int_attribs.resize(finest_level_in_file+1);

    for (int i = 0; i < cnt; i++) {
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
            host_real_attribs[lev][ind][icomp].push_back(*rptr);
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

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt,
          template<class> class Allocator>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt, Allocator>
::WriteAsciiFile (const std::string& filename)
{
    BL_PROFILE("ParticleContainer::WriteAsciiFile()");
    AMREX_ASSERT(!filename.empty());

    const auto strttime = amrex::second();
    //
    // Count # of valid particles.
    //
    Long nparticles = 0;

    for (int lev = 0; lev < m_particles.size();  lev++) {
        auto& pmap = m_particles[lev];
        for (const auto& kv : pmap) {
            const auto& aos = kv.second.GetArrayOfStructs();
            auto np = aos.numParticles();
            Gpu::HostVector<ParticleType> host_aos(np);
            Gpu::copy(Gpu::deviceToHost, aos.begin(), aos.begin() + np, host_aos.begin());
            for (int k = 0; k < np; ++k) {
                const ParticleType& p = host_aos[k];
                if (p.id() > 0)
                    //
                    // Only count (and checkpoint) valid particles.
                    //
                    nparticles++;
            }
        }
    }

    //
    // And send count to I/O processor.
    //
    ParallelDescriptor::ReduceLongSum(nparticles,ParallelDescriptor::IOProcessorNumber());

    if (ParallelDescriptor::IOProcessor())
    {
        //
        // Have I/O processor open file and write out particle metadata.
        //
        std::ofstream File;

        File.open(filename.c_str(), std::ios::out|std::ios::trunc);

        if (!File.good())
            amrex::FileOpenFailed(filename);

        File << nparticles  << '\n';
        File << NStructReal << '\n';
        File << NStructInt  << '\n';
        File << NumRealComps()  << '\n';
        File << NumIntComps()   << '\n';

        File.flush();

        File.close();

        if (!File.good())
            amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");
    }

    ParallelDescriptor::Barrier();

    const int MyProc = ParallelDescriptor::MyProc();

    for (int proc = 0; proc < ParallelDescriptor::NProcs(); proc++)
    {
        if (MyProc == proc)
        {
            //
            // Each CPU opens the file for appending and adds its particles.
            //
            VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

            std::ofstream File;

            File.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

            File.open(filename.c_str(), std::ios::out|std::ios::app);

            File.precision(15);

            if (!File.good())
                amrex::FileOpenFailed(filename);

            for (int lev = 0; lev < m_particles.size();  lev++) {
              auto& pmap = m_particles[lev];
              for (const auto& kv : pmap) {
                ParticleTile<NStructReal, NStructInt, NArrayReal, NArrayInt,
                             amrex::PinnedArenaAllocator> pinned_ptile;
                pinned_ptile.define(NumRuntimeRealComps(), NumRuntimeIntComps());
                pinned_ptile.resize(kv.second.numParticles());
                amrex::copyParticles(pinned_ptile, kv.second);
                const auto& host_aos = pinned_ptile.GetArrayOfStructs();
                const auto& host_soa = pinned_ptile.GetStructOfArrays();

                auto np = host_aos.numParticles();
                for (int index = 0; index < np; ++index) {
                    const ParticleType* it = &host_aos[index];
                    if (it->id() > 0) {

                        // write out the particle struct first...
                        AMREX_D_TERM(File << it->pos(0) << ' ',
                                          << it->pos(1) << ' ',
                                          << it->pos(2) << ' ');

                        for (int i = 0; i < NStructReal; i++)
                            File << it->rdata(i) << ' ';

                        File << it->id()  << ' ';
                        File << it->cpu() << ' ';

                        for (int i = 0; i < NStructInt; i++)
                            File << it->idata(i) << ' ';

                        // then the particle attributes.
                        for (int i = 0; i < NumRealComps(); i++)
                            File << host_soa.GetRealData(i)[index] << ' ';

                        for (int i = 0; i < NumIntComps(); i++)
                            File << host_soa.GetIntData(i)[index] << ' ';

                        File << '\n';
                    }
                }
              }
            }

            File.flush();

            File.close();

            if (!File.good())
                amrex::Abort("ParticleContainer::WriteAsciiFile(): problem writing file");

        }

        ParallelDescriptor::Barrier();
    }

    if (m_verbose > 1)
    {
        auto stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "ParticleContainer::WriteAsciiFile() time: " << stoptime << '\n';
    }
}

#endif /*AMREX_PARTICLEIO_H*/
