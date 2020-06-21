#ifndef AMREX_PARTICLEASYNCIO_H_
#define AMREX_PARTICLEASYNCIO_H_

#include <AMReX_TypeTraits.H>
#include <AMReX_Particles.H>

namespace amrex
{

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::WriteBinaryParticleDataAsync (const std::string& dir, const std::string& name,
                                const Vector<int>& write_real_comp,
                                const Vector<int>& write_int_comp,
                                const Vector<std::string>& real_comp_names,
                                const Vector<std::string>& int_comp_names) const
{
    BL_PROFILE("writeBinaryParticleDataAsync");

    // AMREX_ASSERT(OK());

    // AMREX_ASSERT(sizeof(typename ParticleType::RealType) == 4 ||
    //              sizeof(typename ParticleType::RealType) == 8);

    // const int NProcs = ParallelDescriptor::NProcs();
    // const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

    // AMREX_ALWAYS_ASSERT(real_comp_names.size() == NumRealComps() + NStructReal);
    // AMREX_ALWAYS_ASSERT( int_comp_names.size() == NumIntComps() + NStructInt);

    // std::string pdir = dir;
    // if ( not pdir.empty() and pdir[pdir.size()-1] != '/') pdir += '/';
    // pdir += name;

    // if ( ! levelDirectoriesCreated)
    // {
    //     if (ParallelDescriptor::IOProcessor())
    //         if ( ! amrex::UtilCreateDirectory(pdir, 0755))
    //             amrex::CreateDirectoryFailed(pdir);
    //     ParallelDescriptor::Barrier();
    // }

    // Vector<Long> np_per_grid;
    // Vector<int> lev_offsets(1, 0);
    // for (int lev = 0; lev <= finestLevel(); lev++)
    // {
    //     lev_offsets.push_back(np_per_grid.size());
    //     auto np_pg_lev = NumberOfParticlesInGrid(lev);
    //     np_per_grid.insert(np_per_grid.end(), np_per_grid_lev.begin(), np_per_grid_lev.end());
    // }

    // std::ofstream HdrFile;

    // Long nparticles = 0;
    // int maxnextid;

    // if(usePrePost)
    // {
    //     nparticles = nparticlesPrePost;
    //     maxnextid  = maxnextidPrePost;
    // }
    // else
    // {
    //     nparticles = 0;
    //     maxnextid  = ParticleType::NextID();

    //     for (int lev = 0; lev < m_particles.size();  lev++)
    //     {
    //         const auto& pmap = m_particles[lev];
    //         for (const auto& kv : pmap)
    //         {
    //             const auto& pflags = particle_io_flags[lev][kv.first];
    //             for (int k = 0; k < kv.second.numParticles(); ++k)
    //             {
    //                 if (pflags[k]) nparticles++;
    //             }
    //         }
    //     }

    //     ParallelDescriptor::ReduceLongSum(nparticles, IOProcNumber);
    //     ParticleType::NextID(maxnextid);
    //     ParallelDescriptor::ReduceIntMax(maxnextid, IOProcNumber);
    // }

    // if (ParallelDescriptor::IOProcessor())
    // {
    //     std::string HdrFileName = pdir;

    //     if ( ! HdrFileName.empty() && HdrFileName[HdrFileName.size()-1] != '/')
    //         HdrFileName += '/';

    //     HdrFileName += "Header";
    //     HdrFileNamePrePost = HdrFileName;

    //     HdrFile.open(HdrFileName.c_str(), std::ios::out|std::ios::trunc);

    //     if ( ! HdrFile.good()) amrex::FileOpenFailed(HdrFileName);

    //     //
    //     // First thing written is our Checkpoint/Restart version string.
    //     // We append "_single" or "_double" to the version string indicating
    //     // whether we're using "float" or "double" floating point data in the
    //     // particles so that we can Restart from the checkpoint files.
    //     //
    //     if (sizeof(typename ParticleType::RealType) == 4)
    //     {
    //         HdrFile << ParticleType::Version() << "_single" << '\n';
    //     }
    //     else
    //     {
    //         HdrFile << ParticleType::Version() << "_double" << '\n';
    //     }

    //     int num_output_real = 0;
    //     for (int i = 0; i < NumRealComps() + NStructReal; ++i)
    //         if (write_real_comp[i]) ++num_output_real;

    //     int num_output_int = 0;
    //     for (int i = 0; i < NumIntComps() + NStructInt; ++i)
    //         if (write_int_comp[i]) ++num_output_int;

    //     // AMREX_SPACEDIM and N for sanity checking.
    //     HdrFile << AMREX_SPACEDIM << '\n';

    //     // The number of extra real parameters
    //     HdrFile << num_output_real << '\n';

    //     // Real component names
    //     for (int i = 0; i < NStructReal + NumRealComps(); ++i )
    //         if (write_real_comp[i]) HdrFile << real_comp_names[i] << '\n';

    //     // The number of extra int parameters
    //     HdrFile << num_output_int << '\n';

    //     // int component names
    //     for (int i = 0; i < NStructInt + NumIntComps(); ++i )
    //         if (write_int_comp[i]) HdrFile << int_comp_names[i] << '\n';

    //     bool is_checkpoint = true; // legacy
    //     HdrFile << is_checkpoint << '\n';

    //     // The total number of particles.
    //     HdrFile << nparticles << '\n';

    //     // The value of nextid that we need to restore on restart.
    //     HdrFile << maxnextid << '\n';

    //     // Then the finest level of the AMR hierarchy.
    //     HdrFile << finestLevel() << '\n';

    //     // Then the number of grids at each level.
    //     for (int lev = 0; lev <= finestLevel(); lev++)
    //         HdrFile << ParticleBoxArray(lev).size() << '\n';
    // }

    // // We want to write the data out in parallel.
    // // We'll allow up to nOutFiles active writers at a time.
    // int nOutFiles(256);

    // ParmParse pp("particles");
    // pp.query("particles_nfiles",nOutFiles);
    // if(nOutFiles == -1) nOutFiles = NProcs;
    // nOutFiles = std::max(1, std::min(nOutFiles,NProcs));
    // nOutFilesPrePost = nOutFiles;

    // for (int lev = 0; lev <= finestLevel(); lev++)
    // {
    //     bool gotsome;
    //     if(usePrePost)
    //     {
    //         gotsome = (nParticlesAtLevelPrePost[lev] > 0);
    //     }
    //     else
    //     {
    //         gotsome = (NumberOfParticlesAtLevel(lev) > 0);
    //     }

    //     // We store the particles at each level in their own subdirectory.
    //     std::string LevelDir = pdir;

    //     if (gotsome)
    //     {
    //         if ( ! LevelDir.empty() && LevelDir[LevelDir.size()-1] != '/') LevelDir += '/';

    //         LevelDir = amrex::Concatenate(LevelDir + "Level_", lev, 1);

    //         if ( ! levelDirectoriesCreated) {
    //             if (ParallelDescriptor::IOProcessor())
    //                 if ( ! amrex::UtilCreateDirectory(LevelDir, 0755))
    //                     amrex::CreateDirectoryFailed(LevelDir);
    //             //
    //             // Force other processors to wait until directory is built.
    //             //
    //             ParallelDescriptor::Barrier();
    //         }
    //     }

    //     // Write out the header for each particle
    //     if (gotsome and ParallelDescriptor::IOProcessor()) {
    //         std::string HeaderFileName = LevelDir;
    //         HeaderFileName += "/Particle_H";
    //         std::ofstream ParticleHeader(HeaderFileName);

    //         ParticleBoxArray(lev).writeOn(ParticleHeader);
    //         ParticleHeader << '\n';

    //         ParticleHeader.flush();
    //         ParticleHeader.close();
    //     }

    //     MFInfo info;
    //     info.SetAlloc(false);
    //     MultiFab state(ParticleBoxArray(lev),
    //                    ParticleDistributionMap(lev),
    //                    1,0,info);

    //     // We eventually want to write out the file name and the offset
    //     // into that file into which each grid of particles is written.
    //     Vector<int>  which(state.size(),0);
    //     Vector<int > count(state.size(),0);
    //     Vector<Long> where(state.size(),0);

    //     std::string filePrefix(LevelDir);
    //     filePrefix += '/';
    //     filePrefix += ParticleType::DataPrefix();
    //     if(usePrePost) {
    //         filePrefixPrePost[lev] = filePrefix;
    //     }
    //     bool groupSets(false), setBuf(true);

    //     if (gotsome)
    //     {
    //         for(NFilesIter nfi(nOutFiles, filePrefix, groupSets, setBuf); nfi.ReadyToWrite(); ++nfi)
    //         {
    //             std::ofstream& myStream = (std::ofstream&) nfi.Stream();
    //             WriteParticles(lev, myStream, nfi.FileNumber(), which, count, where,
    //                            write_real_comp, write_int_comp, particle_io_flags);
    //         }

    //         if(usePrePost) {
    //             whichPrePost[lev] = which;
    //             countPrePost[lev] = count;
    //             wherePrePost[lev] = where;
    //         } else {
    //             ParallelDescriptor::ReduceIntSum (which.dataPtr(), which.size(), IOProcNumber);
    //             ParallelDescriptor::ReduceIntSum (count.dataPtr(), count.size(), IOProcNumber);
    //             ParallelDescriptor::ReduceLongSum(where.dataPtr(), where.size(), IOProcNumber);
    //         }
    //     }

    //     if (ParallelDescriptor::IOProcessor())
    //     {
    //         if(usePrePost) {
    //             // ---- write to the header and unlink in CheckpointPost
    //         } else {
    //             for (int j = 0; j < state.size(); j++)
    //             {
    //                 HdrFile << which[j] << ' ' << count[j] << ' ' << where[j] << '\n';
    //             }

    //             if (gotsome && doUnlink)
    //             {
    //                 // Unlink any zero-length data files.
    //                 Vector<Long> cnt(nOutFiles,0);

    //                 for (int i = 0, N=count.size(); i < N; i++) {
    //                     cnt[which[i]] += count[i];
    //                 }

    //                 for (int i = 0, N=cnt.size(); i < N; i++)
    //                 {
    //                     if (cnt[i] == 0)
    //                     {
    //                         std::string FullFileName = NFilesIter::FileName(i, filePrefix);
    //                         FileSystem::Remove(FullFileName);
    //                     }
    //                 }
    //             }
    //         }
    //     }
    // }

    // if (ParallelDescriptor::IOProcessor())
    // {
    //     HdrFile.flush();
    //     HdrFile.close();
    //     if ( ! HdrFile.good())
    //     {
    //         amrex::Abort("ParticleContainer::Checkpoint(): problem writing HdrFile");
    //     }
    // }
}

}

#endif
