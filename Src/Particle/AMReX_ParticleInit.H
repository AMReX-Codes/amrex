
template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitFromAsciiFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= NStructReal);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const Real strttime = amrex::second();
    //
    // Number of processes that read from the file.  //
    int NReaders = ParticleType::MaxReaders();
    //
    // Number of chunks we break the redistribution loop into
    // This default is just for the case of NProcs = NReaders = 1.
    //
    int NRedist = 1;

    if (NProcs <= 1024)
    {
       if (NReaders > 1)
          NRedist = 2;
    }
    else if (NProcs <= 4096)
    {
       NReaders = std::max(NReaders,128);
       NRedist = 4;
    }
    else if (NProcs <= 8192)
    {
       NReaders = std::max(NReaders,384);
       NRedist = 32;
    }
    else if (NProcs <= 16384)
    {
       NReaders = std::max(NReaders,512);
       NRedist = 48;
    }

    resizeData();

    IntVect lNrep(AMREX_D_DECL(1,1,1));

    if (Nrep != 0)
        lNrep = *Nrep;

    long how_many      = 0;
    long how_many_read = 0;
    //
    // NReaders will read particles into a AoS.
    // Later they'll each put their chunk into m_particles and call Redistribute().
    //
    Gpu::HostVector<ParticleType> nparticles;

    if (MyProc < NReaders)
    {
        //
        // Only the first NReaders MPI processes read particles.
        //

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);

        std::ifstream ifs;

        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(file.c_str(), std::ios::in);

        if (!ifs.good())
            amrex::FileOpenFailed(file);

        int cnt = 0;

        ifs >> cnt >> std::ws;

	ParticleLocData pld;
        ParticleType p, p_rep;

        const int Chunk = cnt / NReaders;

        for (int i = 0; i < MyProc*Chunk; i++)
        {
            ifs.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            ifs >> std::ws;  // Eat newline.
        }

        if (!ifs.good())
        {
            std::string msg("ParticleContainer::InitFromAsciiFile(");
            msg += file;
            msg += ") failed @ 1";
            amrex::Error(msg.c_str());
        }

        int MyCnt = Chunk;

        if (MyProc == (NReaders - 1))
            //
            // We'll take the remainder.
            //
            MyCnt += cnt % NReaders;

        const Geometry& geom = Geom(0);

        for (int i = 0; i < MyCnt; i++)
        {
            //
            // We don't read in idata.id or idata.cpu.  We'll set those later
            // in a manner to guarantee the global uniqueness of the pair.
            //
            AMREX_D_TERM(ifs >> p.m_rdata.pos[0];,
                   ifs >> p.m_rdata.pos[1];,
                   ifs >> p.m_rdata.pos[2];);

            for (int n = 0; n < extradata; n++)
            {
                ifs >> p.m_rdata.arr[AMREX_SPACEDIM+n];
            }

            if (!ifs.good())
            {
                std::string msg("ParticleContainer::InitFromAsciiFile(");
                msg += file; msg += ") failed @ 2";
                amrex::Error(msg.c_str());
            }

            if (!Where(p, pld))
            {
                PeriodicShift(p);

                if (!Where(p, pld))
                {
                    if (m_verbose) {
                        amrex::AllPrint() << "BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n'
                                          << "BAD PARTICLE POS " 
                                          << AMREX_D_TERM(   p.m_rdata.pos[0],
                                                             << p.m_rdata.pos[1],
                                                             << p.m_rdata.pos[2])
                                          << "\n";
                    }
                    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromAsciiFile(): invalid particle");
                }
            }

            p.m_idata.id  = ParticleType::NextID();
            p.m_idata.cpu = MyProc;

            nparticles.push_back(p);

            how_many++;
            how_many_read++;

            const Real DomSize[AMREX_SPACEDIM] = { AMREX_D_DECL((geom.ProbHi(0)-geom.ProbLo(0))/lNrep[0],
                                                       (geom.ProbHi(1)-geom.ProbLo(1))/lNrep[1],
                                                       (geom.ProbHi(2)-geom.ProbLo(2))/lNrep[2]) };
            int rep[AMREX_SPACEDIM];

#if AMREX_SPACEDIM==3
            for (rep[2] = 1; rep[2] <= lNrep[2]; rep[2]++)
            {
#endif
                for (rep[1] = 1; rep[1] <= lNrep[1]; rep[1]++)
                {
                    for (rep[0] = 1; rep[0] <= lNrep[0]; rep[0]++) 
                    {
                        if (!(AMREX_D_TERM( (rep[0] == 1), && (rep[1] == 1), && (rep[2] == 1) ) ) )
                        {
                            //
                            // Shift the position.
                            //
                            for (int d=0; d<AMREX_SPACEDIM; ++d)
                            {
                                p_rep.m_rdata.pos[d] = p.m_rdata.pos[d] + (rep[d]-1)*DomSize[d];
                            }
                            //
                            // Copy the mass and velocity.
                            //
                            for (int n = 0; n < extradata; n++)
                                p_rep.m_rdata.arr[AMREX_SPACEDIM+n] = p.m_rdata.arr[AMREX_SPACEDIM+n];

                            if (!Where(p_rep, pld))
                            {
                                PeriodicShift(p_rep);
                                if (!Where(p_rep, pld))
                                {
                                    if (m_verbose) {
                                        amrex::AllPrint() << "BAD REPLICATED PARTICLE ID WOULD BE " << ParticleType::NextID() << "\n";
                                    }
                                    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromAsciiFile(): invalid replicated particle");
                                }
                            }
                            //
                            // Increment the particle ID.
                            //
                            p_rep.m_idata.id  = ParticleType::NextID();
                            //
                            // Assign to the same processor for now.
                            //
                            p_rep.m_idata.cpu = MyProc;
   
                            nparticles.push_back(p_rep);
   
                            how_many++;
                        }
                    }
                }
#if AMREX_SPACEDIM==3
            }
#endif
        }
    }
    //
    // We've read in all the particles.
    // Now Redistribute() each chunk separately to minimize memory bloat.
    //
    int NRedist_chunk = NReaders / NRedist;

    ParticleLocData pld;

    for (int nr = 0; nr < NRedist; nr++)
    {
        if (m_verbose > 0) {
            amrex::Print() << "Redistributing from processor " << nr*NRedist_chunk << " to " 
                                                           << (nr+1)*NRedist_chunk-1 << '\n';
        }
        for (int which = nr*NRedist_chunk; which < (nr+1)*NRedist_chunk; which++)
        {
            if (which == MyProc)
            {
                while (!nparticles.empty())
                {
                    ParticleType& p = nparticles.back();
		    Where(p, pld);

                    m_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
    
                    nparticles.pop_back();
                }
            }
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute();
    }
    //
    // Take care of any leftover chunk
    //
    if (m_verbose > 0)
    {
        if (NRedist*NRedist_chunk < NReaders) {
            amrex::Print() << "Redistributing from processor " << NRedist*NRedist_chunk << " to " 
                           << NReaders << '\n';
        }
    }
    for (int which = NRedist*NRedist_chunk; which < NReaders; which++)
    {
        if (which == MyProc)
        {
	    while (!nparticles.empty())
            {
                ParticleType& p = nparticles.back();
		Where(p, pld);

                m_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);

                nparticles.pop_back();
            }
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute();
    }
    //
    // Add up all the particles read in on each processor to get the total number of particles.
    //
    if (m_verbose > 0)
    {
        const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

        long num_particles = how_many; 

        ParallelDescriptor::ReduceLongSum(num_particles, IOProcNumber);

        if (AMREX_D_TERM(lNrep[0] == 1, && lNrep[1] == 1, && lNrep[2] == 1))
        {
            amrex::Print() << "Total number of particles: " << num_particles << '\n';
        }
        else
        {
            long num_particles_read = how_many_read; 

            ParallelDescriptor::ReduceLongSum(num_particles_read, IOProcNumber);

            amrex::Print() << "Replication the domain with vector           "
                           << AMREX_D_TERM(lNrep[0] << " ", << lNrep[1] << " ", << lNrep[2]) << "\n"
                           << "Total number of particles read in          : " << num_particles_read << '\n'
                           << "Total number of particles after replication: " << num_particles      << '\n';
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "InitFromAsciiFile() time: " << runtime << '\n';
    }
}

//
// The format of a binary particle init file:
//
// NP -- The number of particles in the file.   A "long".
// DM -- Our dimension.  Either 1, 2, or 3.     A "int".
// NX -- The amount of "extra" data.            A "int".
// NP*(DM+NX) native floating-point numbers.    A "float" or "double".
//
// Note that there is nothing separating all these values.
// They're packed into the binary file like sardines.
//

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile (const std::string& file,
                                                                                       int                extradata)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitFromBinaryFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= NStructReal);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = amrex::second();
    //
    // The number of MPI processes that read from the file.
    // We limit this to a rather small number since there's a limit
    // to the number of independent I/O channels on most filesystems.
    //
    const int NReaders = ParticleType::MaxReaders();

    BL_ASSERT(NReaders <= NProcs);
    //
    // How many particles each NReaders reads before redistributing.
    //
    const long NPartPerRedist = ParticleType::MaxParticlesPerRead();

    if (m_verbose > 0)
    {
        amrex::Print() << "Reading with " << NReaders << " readers\n"
                       << "Redistributing after every " << NPartPerRedist << " particles for each reader\n";
    }
    //
    // tmp_particles should mirror how m_particles is built.
    // At the end of this routine it'll be the new m_particles.
    //
    Vector<ParticleLevel> tmp_particles;

    tmp_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
    tmp_particles.resize(finestLevel()+1);
    
    resizeData();

    //
    // All the processors need to participate in Redistribute() though.
    //
    int  NX = 0;
    long NP = 0;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    
    std::ifstream ifs;

    //
    // The "set" of MPI processor numbers of the readers.  We want
    // to be able to easily check whether or not we're a reader.
    //
    std::set<int> readers;
    //
    // The same set but in ascending order.
    //
    Vector<int> rprocs(NReaders);

    if (NReaders == NProcs)
    {
        //
        // Just set'm.
        //
        for (int i = 0; i < NProcs; i++)
            rprocs[i] = i;
    }
    else
    {
        //
        // The I/O Proc builds a set of NReader integers in the range: [0,NProcs-1].
        //
        // It then broadcast'm to all MPI procs.
        //
        // We want these to be as evenly distributed over the full set of
        // [0,NProcs-1] MPI processors as possible, so that when reading we
        // minimize the number of readers per Node, and hence can use more
        // of the available Node memory for reading.
        //
        if (ParallelDescriptor::IOProcessor())
        {
            do
            {
                int n = int(amrex::Random() * (NProcs-1));

                BL_ASSERT(n >= 0);
                BL_ASSERT(n < NProcs);

                readers.insert(n);
            }
            while (readers.size() < NReaders);

            BL_ASSERT(readers.size() == rprocs.size());

            int i = 0;

            for (std::set<int>::const_iterator it = readers.begin(), End = readers.end();
                 it != End;
                 ++it, ++i)
            {
                rprocs[i] = *it;
            }
        }

        ParallelDescriptor::Bcast(rprocs.dataPtr(), rprocs.size(), IOProc);
    }

    if (readers.empty())
    {
        //
        // Set readers for non I/O procs.
        //
        readers.insert(rprocs.begin(), rprocs.end());

        BL_ASSERT(readers.size() == rprocs.size());

        BL_ASSERT(rprocs.size() == NReaders);
    }

    int RealSizeInFile = 0;

    if (readers.find(MyProc) != readers.end())
    {
        int DM = 0;

        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(file.c_str(), std::ios::in|std::ios::binary);

        if (!ifs.good())
            amrex::FileOpenFailed(file);

        ifs.read((char*)&NP, sizeof(NP));
        ifs.read((char*)&DM, sizeof(DM));
        ifs.read((char*)&NX, sizeof(NX));
        //
        // NP MUST be positive!
        //
        if (NP <= 0)
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile(): NP <= 0");
        //
        // DM must equal AMREX_SPACEDIM.
        //
        if (DM != AMREX_SPACEDIM)
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile(): DM != AMREX_SPACEDIM");
        //
        // NX MUST be in [0,N].
        //
        if (NX < 0 || NX > NStructReal)
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile(): NX < 0 || NX > N");
        //
        // Can't ask for more data than exists in the file!
        //
        if (extradata > NX)
            amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile(): extradata > NX");
        //
        // Figure out whether we're dealing with floats or doubles.
        //
        // First get our current position.
        //
        const std::streamoff CURPOS = ifs.tellg();
        //
        // Seek to end of file.
        //
        ifs.seekg(0,std::ios::end);
        //
        // ENDPOS - CURPOS should bracket the particle data.
        //
        const std::streamoff ENDPOS = ifs.tellg();

        RealSizeInFile = (ENDPOS - CURPOS) / (NP*(DM+NX));

        BL_ASSERT(RealSizeInFile == sizeof(float) || RealSizeInFile == sizeof(double));
        //
        // Now set stream back to earlier position.
        //
        ifs.seekg(CURPOS, std::ios::beg);
        //
        // Skip to our place in the file.
        //
        int id = 0;
        for ( ; id < NReaders; id++)
            if (rprocs[id] == MyProc)
                break;

        BL_ASSERT(id >= 0 && id < NReaders);

        const std::streamoff NSKIP = id * (NP/NReaders) * (DM+NX) * RealSizeInFile;

        if (NSKIP > 0)
        {
            ifs.seekg(NSKIP, std::ios::cur);
        }

        if (!ifs.good())
        {
            std::string msg("ParticleContainer::InitFromBinaryFile(");
            msg += file;
            msg += ") failed @ 1";
            amrex::Error(msg.c_str());
        }
    }
    //
    // Everyone needs to know NP -- the number of particles in the file.
    //
    ParallelDescriptor::ReduceLongMax(NP);
    //
    // How many particles each reader gets to read.
    //
    long MyCnt = NP / NReaders;

    if (MyProc == rprocs[0])
        //
        // Give any remainder to the first reader.
        //
        MyCnt += NP % NReaders;

    long how_many_redists = NP / (NPartPerRedist*NReaders), how_many_read = 0;

    if (NP % (NPartPerRedist*NReaders)) how_many_redists++;

    Vector<float>  fxtra, fignore;
    Vector<double> dxtra, dignore;

    if (extradata > 0)
    {
        fxtra.resize(extradata);
        dxtra.resize(extradata);
    }

    if ((NX-extradata) > 0)
    {
        fignore.resize(NX-extradata);
        dignore.resize(NX-extradata);
    }

    ParticleLocData pld;

    for (int j = 0; j < how_many_redists; j++)
    {

        Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
        host_particles.reserve(15);
        host_particles.resize(finestLevel()+1);
        
        if (readers.find(MyProc) != readers.end())
        {
            ParticleType p;

            BL_ASSERT(MyCnt > how_many_read);

            const long NRead = std::min((MyCnt-how_many_read), NPartPerRedist);

            for (long i = 0; i < NRead; i++)
            {
                //
                // We don't read in idata.id or idata.cpu.  We'll set those later
                // in a manner to guarantee the global uniqueness of the pair.
                //
                if (RealSizeInFile == sizeof(float))
                {
                    float fpos[AMREX_SPACEDIM];

                    ifs.read((char*)&fpos[0], AMREX_SPACEDIM*sizeof(float));

                    AMREX_D_TERM(p.m_rdata.pos[0] = fpos[0];,
                           p.m_rdata.pos[1] = fpos[1];,
                           p.m_rdata.pos[2] = fpos[2];);

                }
                else if (RealSizeInFile == sizeof(double))
                {
                    double dpos[AMREX_SPACEDIM];

                    ifs.read((char*)&dpos[0], AMREX_SPACEDIM*sizeof(double));

                    AMREX_D_TERM(p.m_rdata.pos[0] = dpos[0];,
                           p.m_rdata.pos[1] = dpos[1];,
                           p.m_rdata.pos[2] = dpos[2];);
                }

                //
                // Read in any "extradata".
                //
                if (extradata > 0)
                {
                    if (RealSizeInFile == sizeof(float))
                    {
                        ifs.read((char*)&fxtra[0], extradata*sizeof(float));

                        for (int ii = 0; ii < extradata; ii++)
                            p.m_rdata.arr[AMREX_SPACEDIM+ii] = fxtra[ii];
                    }
                    else if (RealSizeInFile == sizeof(double))
                    {
                        ifs.read((char*)&dxtra[0], extradata*sizeof(double));

                        for (int ii = 0; ii < extradata; ii++)
                            p.m_rdata.arr[AMREX_SPACEDIM+ii] = dxtra[ii];
                    }
                }
                //
                // Read any remaining data for this particle.
                //
                if ((NX-extradata) > 0)
                {
                    if (RealSizeInFile == sizeof(float))
                    {
                        ifs.read((char*)&fignore[0], (NX-extradata)*sizeof(float));
                    }
                    else if (RealSizeInFile == sizeof(double))
                    {
                        ifs.read((char*)&dignore[0], (NX-extradata)*sizeof(double));
                    }
                }

                if (!ifs.good())
                {
                    std::string msg("ParticleContainer::InitFromBinaryFile(");
                    msg += file;
                    msg += ") failed @ 2";
                    amrex::Error(msg.c_str());
                }

                if (!Where(p, pld))
                {
                    PeriodicShift(p);

                    if (!Where(p, pld))
                    {
                        if (m_verbose) {
                            amrex::AllPrint() << "BAD PARTICLE ID WOULD BE " << ParticleType::NextID() << '\n'
                                              << "BAD PARTICLE POS " 
                                              << AMREX_D_TERM(   p.m_rdata.pos[0],
                                                                 << p.m_rdata.pos[2],
                                                                 << p.m_rdata.pos[3])
                                              << "\n";
                        }
                        amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryFile(): invalid particle");
                    }
                }

                p.m_idata.id  = ParticleType::NextID();
                p.m_idata.cpu = MyProc;

                host_particles[pld.m_lev][std::make_pair(pld.m_grid, pld.m_tile)].push_back(p);
            }

            how_many_read += NRead;
        }

        for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
        {
            for (auto& kv : host_particles[host_lev]) {
                auto grid = kv.first.first;
                auto tile = kv.first.second;
                const auto& src_tile = kv.second;
                
                auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
                auto old_size = dst_tile.GetArrayOfStructs().size();
                auto new_size = old_size + src_tile.size();
                dst_tile.resize(new_size);
                
                Cuda::thrust_copy(src_tile.begin(),
                                  src_tile.end(),
                                  dst_tile.GetArrayOfStructs().begin() + old_size);
            }
        }
        
        Redistribute();
        //
        // Move particles in m_particles into tmp_particles so that
        // we don't keep trying to redistribute particles that have
        // already been redistributed correctly.
        //
        for (int lev = 0; lev < m_particles.size(); lev++)
        {
            auto& pmap     = m_particles[lev];
            auto& tmp_pmap = tmp_particles[lev];

            for (auto kv : pmap) {
                auto& aos = kv.second.GetArrayOfStructs()();
                auto& tmp_aos = tmp_pmap[kv.first].GetArrayOfStructs()();

                tmp_aos.insert(tmp_aos.end(), aos.begin(), aos.end());
                ParticleVector().swap(aos);
            }

            ParticleLevel().swap(pmap);
        }
    }
    //
    // Make tmp_particles the new m_particles.
    //
    tmp_particles.swap(m_particles);
    //
    // Add up all the particles read in to get the total number of particles.
    //
    if (m_verbose > 0)
    {
        long num_particles_read = how_many_read; 

        ParallelDescriptor::ReduceLongSum(num_particles_read, ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "\nTotal number of particles: " << num_particles_read << '\n';
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "InitFromBinaryFile() time: " << runtime << '\n';
    }

    Gpu::Device::streamSynchronize();
}

//
// This function expects to read a file containing the pathnames of
// binary particles files needing to be read in for input.  It expects
// one file name per line.
//

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitFromBinaryMetaFile (const std::string& metafile,
                                                       int                extradata)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitFromBinaryMetaFile()");
    const Real strttime = amrex::second();

    std::ifstream ifs(metafile.c_str(), std::ios::in);

    std::string file;

    for (;;)
    {
        std::getline(ifs,file);

        if (!ifs.good()) break;

        if (m_verbose > 1)
            amrex::Print() << "InitFromBinaryMetaFile: processing file: " << file << '\n';

        InitFromBinaryFile(file, extradata);
    }

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        amrex::Print() << "InitFromBinaryMetaFile() time: " << runtime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
InitRandom (long                    icount,
            unsigned long           iseed,
            const ParticleInitData& pdata,
            bool                    serialize,
            RealBox                 containing_bx)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitRandom()");
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_gdb != 0);

    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       NProcs   = ParallelDescriptor::NProcs();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = amrex::second();
    const Geometry& geom     = Geom(0);

    Real r, x, len[AMREX_SPACEDIM] = { AMREX_D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };

    // We will enforce that the particles are within the containing_bx.  
    // If containing_bx is not passed in, it defaults to the full domain.
    if (!containing_bx.ok()) containing_bx = geom.ProbDomain();
    
    // containing_bx is assumed to lie within the domain.
    if (!geom.ProbDomain().contains(containing_bx))
    {
        containing_bx.setLo(geom.ProbLo());
        containing_bx.setHi(geom.ProbHi());
    }

    const Real* xlo = containing_bx.lo();
    const Real* xhi = containing_bx.hi();

    amrex::InitRandom(iseed+MyProc);

    if (serialize)
    {
        //
        // We'll let IOProc generate the particles so we get the same
        // positions no matter how many CPUs we have.  This is here
        // mainly for debugging purposes.  It's not really useful for
        // very large numbers of particles.
        //
        //
        Vector<typename ParticleType::RealType> pos(icount*AMREX_SPACEDIM);

        if (ParallelDescriptor::IOProcessor())
        {
            for (long j = 0; j < icount; j++)
            {
                for (int i = 0; i < AMREX_SPACEDIM; i++)
                {
                    do
                    {
                        r = amrex::Random();
                        x = geom.ProbLo(i) + (r * len[i]);
                    }
                    while (x < xlo[i] || x > xhi[i]);

                    pos[j*AMREX_SPACEDIM + i] = x;
                }
            }
        }

        ParallelDescriptor::Bcast(pos.dataPtr(), icount*AMREX_SPACEDIM, IOProc);

        ParticleLocData pld;

        int cnt = 0;

        Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
        host_particles.reserve(15);
        host_particles.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<Real>, NArrayReal > > > host_real_attribs;
        host_real_attribs.reserve(15);
        host_real_attribs.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<int>, NArrayInt > > > host_int_attribs;
        host_int_attribs.reserve(15);
        host_int_attribs.resize(finestLevel()+1);

        for (long j = 0; j < icount; j++)
        {
            ParticleType p;

            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                p.m_rdata.pos[i] = pos[j*AMREX_SPACEDIM + i];
            }

            for (int i = 0; i < NStructReal; i++) {
                p.m_rdata.arr[AMREX_SPACEDIM + i] = pdata.real_struct_data[i];
            }

            if (!Where(p, pld)) {
                amrex::Abort("ParticleContainer::InitRandom(): invalid particle");
            }

            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= finestLevel());
            std::pair<int, int> ind(pld.m_grid, pld.m_tile); 

            const int who = ParticleDistributionMap(pld.m_lev)[pld.m_grid];

            if (who == MyProc) {

                // We own it. Add it at the appropriate level.
                p.m_idata.id  = ParticleType::NextID();
                p.m_idata.cpu = MyProc;

                for (int i = 0; i < NStructInt; i++) {
                    p.m_idata.arr[2 + i] = pdata.int_struct_data[i];
                }
                
                // add the struct
                host_particles[pld.m_lev][ind].push_back(p);

                // add the real...
                for (int i = 0; i < NArrayReal; i++) {
                    host_real_attribs[pld.m_lev][ind][i].push_back(pdata.real_array_data[i]);
		}

                // ... and int array data
                for (int i = 0; i < NArrayInt; i++) {
                    host_int_attribs[pld.m_lev][ind][i].push_back(pdata.int_array_data[i]);
                }

                cnt++;
	      }
        }

	for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
        {
            for (auto& kv : host_particles[host_lev]) {
                auto grid = kv.first.first;
                auto tile = kv.first.second;
                const auto& src_tile = kv.second;
                
                auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
                auto old_size = dst_tile.GetArrayOfStructs().size();
                auto new_size = old_size + src_tile.size();
                dst_tile.resize(new_size);
                
                Cuda::thrust_copy(src_tile.begin(),
                                  src_tile.end(),
                                  dst_tile.GetArrayOfStructs().begin() + old_size);

		for (int i = 0; i < NArrayReal; ++i) {
		  Cuda::thrust_copy(host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
		}

		for (int i = 0; i < NArrayInt; ++i) {
		  Cuda::thrust_copy(host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
		}
            }
        }

        BL_ASSERT(OK());
    }
    else {
        // We'll generate the particles in parallel.
        // Each CPU will key off the given seed to get independent streams of random numbers.
        long M = icount / NProcs;
        // Processor 0 will get the slop.
        if (MyProc == 0) {
            M += (icount % NProcs);
        }
        
        ParticleLocData pld;

        Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
        host_particles.reserve(15);
        host_particles.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<Real>, NArrayReal > > > host_real_attribs;
        host_real_attribs.reserve(15);
        host_real_attribs.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<int>, NArrayInt > > > host_int_attribs;
        host_int_attribs.reserve(15);
        host_int_attribs.resize(finestLevel()+1);

        for (long icnt = 0; icnt < M; icnt++) {
            ParticleType p;
            for (int i = 0; i < AMREX_SPACEDIM; i++) {
                do {
                    r = amrex::Random();
                    x = geom.ProbLo(i) + (r * len[i]);
                }
                while (x < xlo[i] || x > xhi[i]);
                
                p.m_rdata.pos[i] = x;
                
                BL_ASSERT(p.m_rdata.pos[i] < geom.ProbHi(i));
            }
            
            for (int i = 0; i < NStructReal; i++) {
                p.m_rdata.arr[AMREX_SPACEDIM + i] = pdata.real_struct_data[i];
            }
            
            // the int struct data
            p.m_idata.id  = ParticleType::NextID();
            p.m_idata.cpu = ParallelDescriptor::MyProc();
            
            for (int i = 0; i < NStructInt; i++) {
                p.m_idata.arr[2 + i] = pdata.int_struct_data[i];
            }
            
            // locate the particle
            if (!Where(p, pld))
            {
                amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitRandom(): invalid particle");
            }
            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= finestLevel());            
            std::pair<int, int> ind(pld.m_grid, pld.m_tile); 
            

	    // add the struct
	    host_particles[pld.m_lev][ind].push_back(p);

	    // add the real...
	    for (int i = 0; i < NArrayReal; i++) {
	        host_real_attribs[pld.m_lev][ind][i].push_back(pdata.real_array_data[i]);
	    }

	    // ... and int array data
	    for (int i = 0; i < NArrayInt; i++) {
	        host_int_attribs[pld.m_lev][ind][i].push_back(pdata.int_array_data[i]);
	    }            
        }

        for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
        {
            for (auto& kv : host_particles[host_lev]) {
                auto grid = kv.first.first;
                auto tile = kv.first.second;
                const auto& src_tile = kv.second;
                
                auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
                auto old_size = dst_tile.GetArrayOfStructs().size();
                auto new_size = old_size + src_tile.size();
                dst_tile.resize(new_size);
                
                Cuda::thrust_copy(src_tile.begin(),
                                  src_tile.end(),
                                  dst_tile.GetArrayOfStructs().begin() + old_size);

		for (int i = 0; i < NArrayReal; ++i) {
		  Cuda::thrust_copy(host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
		}

		for (int i = 0; i < NArrayInt; ++i) {
		  Cuda::thrust_copy(host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
		}
            }
        }
        // Let Redistribute() sort out where the particles belong.
        Redistribute();
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitRandom() time: " << stoptime << '\n';
    }

    Gpu::Device::streamSynchronize();
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>
::InitRandomPerBox (long                    icount_per_box,
                    unsigned long           iseed,
                    const ParticleInitData& pdata)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitRandomPerBox()");
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount_per_box > 0);

    BL_ASSERT(m_gdb != 0);

    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = amrex::second();
    const Geometry& geom     = Geom(0);

    ParticleLocData pld;
    ParticleType p;

    // This assumes level 0 since geom = m_gdb->Geom(0)
    const Real* dx  = geom.CellSize();

    // We use exactly the same seed for every grid
    std::mt19937 mt(iseed);
    std::uniform_real_distribution<double> dist(0.0, 1.0);

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    // We'll generate the particles in parallel -- but no tiling here.
    for (MFIter mfi(*m_dummy_mf[0], false); mfi.isValid(); ++mfi)
    {
        Box grid = m_gdb->ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box = RealBox(grid,dx,geom.ProbLo());

        for (long icnt = 0; icnt < icount_per_box; icnt++) {
        for (long jcnt = 0; jcnt < icount_per_box; jcnt++) {
        for (long kcnt = 0; kcnt < icount_per_box; kcnt++)
        {
            p.pos(0) = grid_box.lo(0) + (dist(mt) + icnt) / icount_per_box;
            p.pos(1) = grid_box.lo(1) + (dist(mt) + jcnt) / icount_per_box;
            p.pos(2) = grid_box.lo(2) + (dist(mt) + kcnt) / icount_per_box;
            
            for (int i = 0; i < AMREX_SPACEDIM; i++)
                BL_ASSERT(p.pos(i) < grid_box.hi(i));

            // the real struct data
            for (int i = 0; i < NStructReal; i++) {
                p.m_rdata.arr[AMREX_SPACEDIM + i] = pdata.real_struct_data[i];
            }

            // the int struct data
            p.m_idata.id  = ParticleType::NextID();
            p.m_idata.cpu = ParallelDescriptor::MyProc();
            
            for (int i = 0; i < NStructInt; i++) {
                p.m_idata.arr[2 + i] = pdata.int_struct_data[i];
            }

            // locate the particle
            if (!Where(p, pld)) {
                amrex::Abort("ParticleContainer::InitRandomPerBox(): invalid particle");
            }
            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= finestLevel());
            std::pair<int, int> ind(pld.m_grid, pld.m_tile); 

            // add the struct
            m_particles[pld.m_lev][ind].push_back(p);
            
            // add the real...
            for (int i = 0; i < NArrayReal; i++) {
                m_particles[pld.m_lev][ind].push_back_real(i, pdata.real_array_data[i]);
            }
            
            // ... and int array data
            for (int i = 0; i < NArrayInt; i++) {
                m_particles[pld.m_lev][ind].push_back_int(i, pdata.int_array_data[i]);
            }            

        } } }
    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitRandomPerBox() time: " << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
InitOnePerCell (Real x_off, Real y_off, Real z_off, const ParticleInitData& pdata)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitOnePerCell()");

    BL_ASSERT(m_gdb != 0);

    // Note that x_off, y_off and z_off are the offsets from the lower left corner
    // in each cell as measured in units of dx, so they must be in [0,1]
    BL_ASSERT(x_off >= 0. && y_off >= 0. && z_off >= 0.); 
    BL_ASSERT(x_off <= 1. && y_off <= 1. && z_off <= 1.); 

    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = amrex::second();
    const Geometry& geom     = Geom(0);

    // This assumes level 0 since geom = Geom(0)
    const Real* dx  = geom.CellSize();

    ParticleLocData pld;
    ParticleType p;

    // We'll generate the particles in parallel -- but no tiling of the grid here.
    for (MFIter mfi(*m_dummy_mf[0], false); mfi.isValid(); ++mfi) {
        Box grid = ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box (grid,dx,geom.ProbLo());
        
        for (IntVect beg = grid.smallEnd(), end=grid.bigEnd(),
                 cell = grid.smallEnd(); cell <= end; grid.next(cell)) { 

            // the real struct data
            AMREX_D_TERM(p.m_rdata.pos[0] = grid_box.lo(0) + (x_off + cell[0]-beg[0])*dx[0];,
                   p.m_rdata.pos[1] = grid_box.lo(1) + (y_off + cell[1]-beg[1])*dx[1];,
                   p.m_rdata.pos[2] = grid_box.lo(2) + (z_off + cell[2]-beg[2])*dx[2];);
            
            for (int d = 0; d < AMREX_SPACEDIM; ++d) {
                BL_ASSERT(p.m_rdata.pos[d] < grid_box.hi(d));
            }

            for (int i = 0; i < NStructReal; i++) {
                p.m_rdata.arr[AMREX_SPACEDIM + i] = pdata.real_struct_data[i];
            }

            // the int struct data
            p.m_idata.id  = ParticleType::NextID();
            p.m_idata.cpu = ParallelDescriptor::MyProc();

            for (int i = 0; i < NStructInt; i++) {
                p.m_idata.arr[2 + i] = pdata.int_struct_data[i];
            }

            // locate the particle
            if (!Where(p, pld))
            {
                amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitOnePerCell(): invalid particle");
            }
            BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= finestLevel());            
            std::pair<int, int> ind(pld.m_grid, pld.m_tile); 

            // add the struct
	    m_particles[pld.m_lev][ind].push_back(p);

            // add the real...
            for (int i = 0; i < NArrayReal; i++) {
                m_particles[pld.m_lev][ind].push_back_real(i, pdata.real_array_data[i]);
            }

            // ... and int array data
            for (int i = 0; i < NArrayInt; i++) {
                m_particles[pld.m_lev][ind].push_back_int(i, pdata.int_array_data[i]);
            }            
        }
    }
    
    if (m_verbose > 1) {
        Real stoptime = amrex::second() - strttime;
        
        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitOnePerCell() time: " << stoptime << '\n';
    }
}

template <int NStructReal, int NStructInt, int NArrayReal, int NArrayInt>
void
ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::
InitNRandomPerCell (int n_per_cell, const ParticleInitData& pdata)
{
    BL_PROFILE("ParticleContainer<NSR, NSI, NAR, NAI>::InitNRandomPerCell()");

    BL_ASSERT(m_gdb != 0);

    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = amrex::second();
    const Geometry& geom     = Geom(0);

    // This assumes level 0 since geom = Geom(0)
    const Real* dx  = geom.CellSize();

    resizeData();

    for (int lev = 0; lev < m_particles.size(); lev++)
        BL_ASSERT(m_particles[lev].empty());

    ParticleLocData pld;
    ParticleType p;
    Real r;

    // We'll generate the particles in parallel -- but no tiling here.
    for (MFIter mfi(*m_dummy_mf[0], false); mfi.isValid(); ++mfi)
    {
        Box grid = ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box (grid,dx,geom.ProbLo());

        Vector<std::map<std::pair<int, int>, Gpu::HostVector<ParticleType> > > host_particles;
        host_particles.reserve(15);
        host_particles.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<Real>, NArrayReal > > > host_real_attribs;
        host_real_attribs.reserve(15);
        host_real_attribs.resize(finestLevel()+1);

	Vector<std::map<std::pair<int, int>, std::array<Cuda::HostVector<int>, NArrayInt > > > host_int_attribs;
        host_int_attribs.reserve(15);
        host_int_attribs.resize(finestLevel()+1);
        
        for (IntVect beg = grid.smallEnd(), end=grid.bigEnd(),
                    cell = grid.smallEnd(); cell <= end; grid.next(cell)) { 

            for (int n = 0; n < n_per_cell; n++)
            {
                // the real struct data
                for (int i = 0; i < AMREX_SPACEDIM; i++) {
                    r = amrex::Random();
                    p.m_rdata.pos[i] = grid_box.lo(i) + (r + cell[i]-beg[i])*dx[i];
                    BL_ASSERT(p.m_rdata.pos[i] < grid_box.hi(i));
                }
                                
                for (int i = 0; i < NStructReal; i++) {
                    p.m_rdata.arr[AMREX_SPACEDIM + i] = pdata.real_struct_data[i];
                }

                // the int struct data
                p.m_idata.id  = ParticleType::NextID();
                p.m_idata.cpu = ParallelDescriptor::MyProc();

                for (int i = 0; i < NStructInt; i++) {
                    p.m_idata.arr[2 + i] = pdata.int_struct_data[i];
                }

                // locate the particle
                if (!Where(p, pld)) {
                    amrex::Abort("ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitNRandomPerCell(): invalid particle");
                }
                BL_ASSERT(pld.m_lev >= 0 && pld.m_lev <= finestLevel());
                std::pair<int, int> ind(pld.m_grid, pld.m_tile); 

                // add the struct
                host_particles[pld.m_lev][ind].push_back(p);

                // add the real...
                for (int i = 0; i < NArrayReal; i++) {
                    host_real_attribs[pld.m_lev][ind][i].push_back(pdata.real_array_data[i]);
		}

                // ... and int array data
                for (int i = 0; i < NArrayInt; i++) {
                    host_int_attribs[pld.m_lev][ind][i].push_back(pdata.int_array_data[i]);
                }
            }
        }

	for (int host_lev = 0; host_lev < static_cast<int>(host_particles.size()); ++host_lev)
        {
            for (auto& kv : host_particles[host_lev]) {
                auto grid = kv.first.first;
                auto tile = kv.first.second;
                const auto& src_tile = kv.second;
                
                auto& dst_tile = GetParticles(host_lev)[std::make_pair(grid,tile)];
                auto old_size = dst_tile.GetArrayOfStructs().size();
                auto new_size = old_size + src_tile.size();
                dst_tile.resize(new_size);
                
                Cuda::thrust_copy(src_tile.begin(),
                                  src_tile.end(),
                                  dst_tile.GetArrayOfStructs().begin() + old_size);

		for (int i = 0; i < NArrayReal; ++i) {
		  Cuda::thrust_copy(host_real_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_real_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetRealData(i).begin() + old_size);
		}

		for (int i = 0; i < NArrayInt; ++i) {
		  Cuda::thrust_copy(host_int_attribs[host_lev][std::make_pair(grid,tile)][i].begin(),
				    host_int_attribs[host_lev][std::make_pair(grid,tile)][i].end(),
				    dst_tile.GetStructOfArrays().GetIntData(i).begin() + old_size);
		}
            }
        }

    }

    if (m_verbose > 1)
    {
        Real stoptime = amrex::second() - strttime;
        
        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);
        
        amrex::Print() << "ParticleContainer<NStructReal, NStructInt, NArrayReal, NArrayInt>::InitNRandomPerCell() time: " << stoptime << '\n';
    }

    Gpu::Device::streamSynchronize();
}
