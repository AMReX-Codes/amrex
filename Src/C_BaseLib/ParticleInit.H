template <int N>
void
ParticleContainer<N>::InitFromAsciiFile (const std::string& file, int extradata, const IntVect* Nrep)
{
    BL_PROFILE("ParticleContainer<N>::InitFromAsciiFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= N);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const Real strttime = ParallelDescriptor::second();
    //
    // Number of processes that read from the file.  //
    int NReaders = ParticleBase::MaxReaders();
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

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
    m_particles.resize(m_gdb->finestLevel()+1);

    IntVect lNrep(D_DECL(1,1,1));

    if (Nrep != 0)
        lNrep = *Nrep;

    long how_many      = 0;
    long how_many_read = 0;
    //
    // NReaders will read particles into a PBox.
    // Later they'll each put their chunk into m_particles and call Redistribute().
    //
    PBox nparticles;

    if (MyProc < NReaders)
    {
        //
        // Only the first NReaders MPI processes read particles.
        //
        std::ifstream ifs;

        VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
 
        ifs.rdbuf()->pubsetbuf(io_buffer.dataPtr(), io_buffer.size());

        ifs.open(file.c_str(), std::ios::in);

        if (!ifs.good())
            BoxLib::FileOpenFailed(file);

        int cnt = 0;

        ifs >> cnt >> std::ws;

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
            BoxLib::Error(msg.c_str());
        }

        int MyCnt = Chunk;

        if (MyProc == (NReaders - 1))
            //
            // We'll take the remainder.
            //
            MyCnt += cnt % NReaders;

        const Geometry& geom = m_gdb->Geom(0);

        const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
        const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
        const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                                  Real(.125)*geom.CellSize(1),
                                                  Real(.125)*geom.CellSize(2)) };
        for (int i = 0; i < MyCnt; i++)
        {
            //
            // We don't read in m_id or m_cpu.  We'll set those later
            // in a manner to guarantee the global uniqueness of the pair.
            //
            D_TERM(ifs >> p.m_pos[0];,
                   ifs >> p.m_pos[1];,
                   ifs >> p.m_pos[2];);
            //
            // If particle is right on a domain boundary then move it just inside the boundary.
            //
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                if (p.m_pos[d] <= ProbLo[d]) p.m_pos[d] += Delta[d];
                if (p.m_pos[d] >= ProbHi[d]) p.m_pos[d] -= Delta[d];
            }

            for (int n = 0; n < extradata; n++)
            {
                ifs >> p.m_data[n];
            }

            if (!ifs.good())
            {
                std::string msg("ParticleContainer::InitFromAsciiFile(");
                msg += file; msg += ") failed @ 2";
                BoxLib::Error(msg.c_str());
            }

            if (!ParticleBase::Where(p,m_gdb))
            {
                ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                {
                    std::cout << "BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                    for (int d = 0; d < BL_SPACEDIM; d++)
                    {
                        std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                    }

                    BoxLib::Abort("ParticleContainer<N>::InitFromAsciiFile(): invalid particle");
                }
            }

            p.m_id  = ParticleBase::NextID();
            p.m_cpu = MyProc;

            nparticles.push_back(p);

            how_many++;
            how_many_read++;

            const Real DomSize[BL_SPACEDIM] = { D_DECL((geom.ProbHi(0)-geom.ProbLo(0))/lNrep[0],
                                                       (geom.ProbHi(1)-geom.ProbLo(1))/lNrep[1],
                                                       (geom.ProbHi(2)-geom.ProbLo(2))/lNrep[2]) };
            int rep[BL_SPACEDIM];

#if BL_SPACEDIM==3
            for (rep[2] = 1; rep[2] <= lNrep[2]; rep[2]++)
            {
#endif
                for (rep[1] = 1; rep[1] <= lNrep[1]; rep[1]++)
                {
                    for (rep[0] = 1; rep[0] <= lNrep[0]; rep[0]++) 
                    {
                        if (!(D_TERM( (rep[0] == 1), && (rep[1] == 1), && (rep[2] == 1) ) ) )
                        {
                            //
                            // Shift the position.
                            //
                            for (int d=0; d<BL_SPACEDIM; ++d)
                            {
                                p_rep.m_pos[d] = p.m_pos[d] + (rep[d]-1)*DomSize[d];
                            }
                            //
                            // Copy the mass and velocity.
                            //
                            for (int n = 0; n < extradata; n++)
                                p_rep.m_data[n] = p.m_data[n];

                            if (!ParticleBase::Where(p_rep,m_gdb))
                            {
                                ParticleBase::PeriodicShift(p_rep,m_gdb);
                                if (!ParticleBase::Where(p_rep,m_gdb))
                                {
                                    std::cout << "BAD REPLICATED PARTICLE ID WOULD BE " << ParticleBase::NextID() << std::endl;
                                    BoxLib::Abort("ParticleContainer<N>::InitFromAsciiFile(): invalid replicated particle");
                                }
                            }
                            //
                            // Increment the particle ID.
                            //
                            p_rep.m_id  = ParticleBase::NextID();
                            //
                            // Assign to the same processor for now.
                            //
                            p_rep.m_cpu = MyProc;
   
                            nparticles.push_back(p_rep);
   
                            how_many++;
                        }
                    }
                }
#if BL_SPACEDIM==3
            }
#endif
        }
    }
    //
    // We've read in all the particles.
    // Now Redistribute() each chunk separately to minimize memory bloat.
    //
    int NRedist_chunk = NReaders / NRedist;

    for (int nr = 0; nr < NRedist; nr++)
    {
        if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
            std::cout << "Redistributing from processor " << nr*NRedist_chunk << " to " 
                                                          << (nr+1)*NRedist_chunk-1 << '\n';
        for (int which = nr*NRedist_chunk; which < (nr+1)*NRedist_chunk; which++)
        {
            if (which == MyProc)
            {
                while (!nparticles.empty())
                {
                    ParticleType& p = nparticles.back();
    
                    m_particles[p.m_lev][p.m_grid].push_back(p);
    
                    nparticles.pop_back();
                }

                PBox().swap(nparticles);
            }
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute(true);
    }
    //
    // Take care of any leftover chunk
    //
    if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        if (NRedist*NRedist_chunk < NReaders)
            std::cout << "Redistributing from processor " << NRedist*NRedist_chunk << " to " 
                                                          << NReaders << '\n';
    }
    for (int which = NRedist*NRedist_chunk; which < NReaders; which++)
    {
        if (which == MyProc)
        {
            while (!nparticles.empty())
            {
                ParticleType& p = nparticles.back();

                m_particles[p.m_lev][p.m_grid].push_back(p);

                nparticles.pop_back();
            }

            PBox().swap(nparticles);
        }
        //
        // Let Redistribute() sort'm out.
        //
        Redistribute(true);
    }
    //
    // Add up all the particles read in on each processor to get the total number of particles.
    //
    if (m_verbose > 0)
    {
        const int IOProcNumber = ParallelDescriptor::IOProcessorNumber();

        long num_particles = how_many; 

        ParallelDescriptor::ReduceLongSum(num_particles, IOProcNumber);

        if (D_TERM(lNrep[0] == 1, && lNrep[1] == 1, && lNrep[2] == 1))
        {
            if (ParallelDescriptor::IOProcessor())
                std::cout << "Total number of particles: " << num_particles << '\n';
        }
        else
        {
            long num_particles_read = how_many_read; 

            ParallelDescriptor::ReduceLongSum(num_particles_read, IOProcNumber);

            if (ParallelDescriptor::IOProcessor())
            {
                std::cout << "Replication the domain with vector           ";
                for (int d=0; d<BL_SPACEDIM; ++d) {
                    std::cout << lNrep[d] << " ";
                }
                std::cout << '\n';
                std::cout << "Total number of particles read in          : " << num_particles_read << '\n';
                std::cout << "Total number of particles after replication: " << num_particles      << '\n';
            }
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromAsciiFile() time: " << runtime << '\n';
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

template <int N>
void
ParticleContainer<N>::InitFromBinaryFile (const std::string& file,
                                          int                extradata)
{
    BL_PROFILE("ParticleContainer<N>::InitFromBinaryFile()");
    BL_ASSERT(!file.empty());
    BL_ASSERT(extradata <= N);

    const int  MyProc   = ParallelDescriptor::MyProc();
    const int  NProcs   = ParallelDescriptor::NProcs();
    const int  IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real strttime = ParallelDescriptor::second();
    //
    // The number of MPI processes that read from the file.
    // We limit this to a rather small number since there's a limit
    // to the number of independent I/O channels on most filesystems.
    //
    const int NReaders = ParticleBase::MaxReaders();

    BL_ASSERT(NReaders <= NProcs);
    //
    // How many particles each NReaders reads before redistributing.
    //
    const long NPartPerRedist = ParticleBase::MaxParticlesPerRead();

    if (m_verbose > 0 && ParallelDescriptor::IOProcessor())
    {
        std::cout << "Reading with " << NReaders << " readers\n";
        std::cout << "Redistributing after every " << NPartPerRedist << " particles for each reader\n";
    }
    //
    // tmp_particles should mirror how m_particles is built.
    // At the end of this routine it'll be the new m_particles.
    //
    Array<PMap> tmp_particles;

    tmp_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
      m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    tmp_particles.resize(m_gdb->finestLevel()+1);
      m_particles.resize(m_gdb->finestLevel()+1);
    //
    // All the processors need to participate in Redistribute() though.
    //
    int  NX = 0;
    long NP = 0;

    std::ifstream ifs;

    VisMF::IO_Buffer io_buffer(VisMF::IO_Buffer_Size);
    //
    // The "set" of MPI processor numbers of the readers.  We want
    // to be able to easily check whether or not we're a reader.
    //
    std::set<int> readers;
    //
    // The same set but in ascending order.
    //
    Array<int> rprocs(NReaders);

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
                int n = int(BoxLib::Random() * (NProcs-1));

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
            BoxLib::FileOpenFailed(file);

        ifs.read((char*)&NP, sizeof(NP));
        ifs.read((char*)&DM, sizeof(DM));
        ifs.read((char*)&NX, sizeof(NX));
        //
        // NP MUST be positive!
        //
        if (NP <= 0)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): NP <= 0");
        //
        // DM must equal BL_SPACEDIM.
        //
        if (DM != BL_SPACEDIM)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): DM != BL_SPACEDIM");
        //
        // NX MUST be in [0,N].
        //
        if (NX < 0 || NX > N)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): NX < 0 || NX > N");
        //
        // Can't ask for more data than exists in the file!
        //
        if (extradata > NX)
            BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): extradata > NX");
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
            BoxLib::Error(msg.c_str());
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

    const Geometry& geom = m_gdb->Geom(0);

    const Real ProbLo[BL_SPACEDIM] = { D_DECL(geom.ProbLo(0), geom.ProbLo(1), geom.ProbLo(2)) };
    const Real ProbHi[BL_SPACEDIM] = { D_DECL(geom.ProbHi(0), geom.ProbHi(1), geom.ProbHi(2)) };
    const Real  Delta[BL_SPACEDIM] = { D_DECL(Real(.125)*geom.CellSize(0),
                                              Real(.125)*geom.CellSize(1),
                                              Real(.125)*geom.CellSize(2)) };

    long how_many_redists = NP / (NPartPerRedist*NReaders), how_many_read = 0;

    if (NP % (NPartPerRedist*NReaders)) how_many_redists++;

    Array<float>  fxtra, fignore;
    Array<double> dxtra, dignore;

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

    for (int j = 0; j < how_many_redists; j++)
    {
        if (readers.find(MyProc) != readers.end())
        {
            ParticleType p;

            BL_ASSERT(MyCnt > how_many_read);

            const long NR = std::min((MyCnt-how_many_read), NPartPerRedist);

            for (long i = 0; i < NR; i++)
            {
                //
                // We don't read in m_id or m_cpu.  We'll set those later
                // in a manner to guarantee the global uniqueness of the pair.
                //
                if (RealSizeInFile == sizeof(float))
                {
                    float fpos[BL_SPACEDIM];

                    ifs.read((char*)&fpos[0], BL_SPACEDIM*sizeof(float));

                    D_TERM(p.m_pos[0] = fpos[0];,
                           p.m_pos[1] = fpos[1];,
                           p.m_pos[2] = fpos[2];);

                }
                else if (RealSizeInFile == sizeof(double))
                {
                    double dpos[BL_SPACEDIM];

                    ifs.read((char*)&dpos[0], BL_SPACEDIM*sizeof(double));

                    D_TERM(p.m_pos[0] = dpos[0];,
                           p.m_pos[1] = dpos[1];,
                           p.m_pos[2] = dpos[2];);
                }
                //
                // If particle is right on a domain boundary then move it just inside the boundary.
                //
                for (int d = 0; d < BL_SPACEDIM; d++)
                {
                    if (p.m_pos[d] <= ProbLo[d]) p.m_pos[d] += Delta[d];
                    if (p.m_pos[d] >= ProbHi[d]) p.m_pos[d] -= Delta[d];
                }
                //
                // Read in any "extradata".
                //
                if (extradata > 0)
                {
                    if (RealSizeInFile == sizeof(float))
                    {
                        ifs.read((char*)&fxtra[0], extradata*sizeof(float));

                        for (int i = 0; i < extradata; i++)
                            p.m_data[i] = fxtra[i];
                    }
                    else if (RealSizeInFile == sizeof(double))
                    {
                        ifs.read((char*)&dxtra[0], extradata*sizeof(double));

                        for (int i = 0; i < extradata; i++)
                            p.m_data[i] = dxtra[i];
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
                    BoxLib::Error(msg.c_str());
                }

                if (!ParticleBase::Where(p,m_gdb))
                {
                    ParticleBase::PeriodicShift(p,m_gdb);

                    if (!ParticleBase::Where(p,m_gdb))
                    {
                        std::cout << "BAD PARTICLE ID WOULD BE " << ParticleBase::NextID() << '\n';

                        for (int d = 0; d < BL_SPACEDIM; d++)
                        {
                            std::cout << "BAD PARTICLE POS(" << d << ") " << p.m_pos[d] << std::endl;
                        }

                        BoxLib::Abort("ParticleContainer<N>::InitFromBinaryFile(): invalid particle");
                    }
                }

                p.m_id  = ParticleBase::NextID();
                p.m_cpu = MyProc;

                m_particles[p.m_lev][p.m_grid].push_back(p);
            }

            how_many_read += NR;
        }

        Redistribute(true);
        //
        // Move particles in m_particles into tmp_particles so that
        // we don't keep trying to redistribute particles that have
        // already been redistributed correctly.
        //
        for (int lev = 0; lev < m_particles.size(); lev++)
        {
            PMap& pmap     = m_particles[lev];
            PMap& tmp_pmap = tmp_particles[lev];

            for (typename PMap::iterator pmap_it = pmap.begin(), End = pmap.end();
                 pmap_it != End;
                 ++pmap_it)
            {
                PBox& pbox     = pmap_it->second;
                PBox& tmp_pbox = tmp_pmap[pmap_it->first];

                tmp_pbox.insert(tmp_pbox.end(), pbox.begin(), pbox.end());

                PBox().swap(pbox);
            }

            PMap().swap(pmap);
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

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "\nTotal number of particles: " << num_particles_read << '\n';
        }
    }

    BL_ASSERT(OK());

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryFile() time: " << runtime << '\n';
    }
}

//
// This function expects to read a file containing the pathnames of
// binary particles files needing to be read in for input.  It expects
// one file name per line.
//

template <int N>
void
ParticleContainer<N>::InitFromBinaryMetaFile (const std::string& metafile,
                                              int                extradata)
{
    BL_PROFILE("ParticleContainer<N>::InitFromBinaryMetaFile()");
    const Real strttime = ParallelDescriptor::second();

    std::ifstream ifs(metafile.c_str(), std::ios::in);

    std::string file;

    for (;;)
    {
        std::getline(ifs,file);

        if (!ifs.good()) break;

        if (m_verbose > 1 && ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryMetaFile: processing file: " << file << '\n';

        InitFromBinaryFile(file, extradata);
    }

    if (m_verbose > 1)
    {
        ByteSpread();

        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, ParallelDescriptor::IOProcessorNumber());

        if (ParallelDescriptor::IOProcessor())
            std::cout << "InitFromBinaryMetaFile() time: " << runtime << '\n';
    }
}

template <int N>
void
ParticleContainer<N>::InitCosmo1ppcMultiLevel(MultiFab& mf, const Real disp_fac[], const Real vel_fac[], 
                                              const Real particleMass, int disp_idx, int vel_idx, 
                                              BoxArray &baWhereNot, int lev)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppcMultiLevel()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(lev);
    const Real*     dx       = geom.CellSize();

    static Array<int> calls;

    calls.resize(m_gdb->initialBaLevels()+1);

    calls[lev]++;

    if (calls[lev] > 1) return;

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->initialBaLevels()+1);

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    Real         vel[BL_SPACEDIM];
    
    Real 	mean_disp[BL_SPACEDIM]={D_DECL(0,0,0)};


    //
    // The mf should be initialized according to the ics...
    //
    int outside_counter=0;
    long outcount[3]={0,0,0};
    long outcountminus[3]={0,0,0};
    long totalcount=0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));
		    totalcount++;
		    if (baWhereNot.contains(indices)) 
		    {
                       continue;
		    }

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,disp_idx+n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
			//
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+Real(0.5))*dx[n];
			if(disp[n]*disp_fac[n]>dx[n]/2.0)
			  outcount[n]++;
			if(disp[n]*disp_fac[n]<-dx[n]/2.0)
			  outcountminus[n]++;
			mean_disp[n]+=fabs(disp[n]);
			//
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.m_pos[n] += disp[n] * disp_fac[n];

                        //
			// Set the velocities.
                        //
                        vel[n] = myFab(indices,vel_idx+n);
	                p.m_data[n+1] = vel[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppcMultiLevel():invalid particle");
                    }

		    BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
		    //handle particles that ran out of this level into a finer one. 
		    if (baWhereNot.contains(p.m_cell))
		    {
		      outside_counter++;
		      ParticleType newp[8];
		      for (int i=0;i<8;i++)
		      {
			newp[i].m_data[0] = particleMass/8.0;
			newp[i].m_id = ParticleBase::NextID();
			newp[i].m_cpu     = MyProc;
			for (int dim=0;dim<BL_SPACEDIM;dim++)
			{
			  newp[i].m_pos[dim]=p.m_pos[dim]+(2*((i/(1 << dim)) % 2)-1)*dx[dim]/4.0;
			  newp[i].m_data[dim+1]=p.m_data[dim+1];
			  
			  
			}
		      
			if (!ParticleBase::Where(newp[i],m_gdb))
			{
			    ParticleBase::PeriodicShift(newp[i],m_gdb);
    
			    if (!ParticleBase::Where(newp[i],m_gdb))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppcMultiLevel():invalid particle");
			}
			m_particles[newp[i].m_lev][newp[i].m_grid].push_back(newp[i]);
		      }
		      
		    }
	            
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
		    else
		      m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
    Redistribute();
}

template <int N>
void
ParticleContainer<N>::InitCosmo1ppc (MultiFab& mf, const Real vel_fac[], const Real particleMass)
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo1ppc()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const Geometry& geom     = m_gdb->Geom(0);
    const Real*     dx       = geom.CellSize();

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    ParticleType p;
    Real         disp[BL_SPACEDIM];
    const Real   len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                             geom.ProbLength(1),
                                             geom.ProbLength(2)) };
    //
    // The grid should be initialized according to the ics...
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        FArrayBox&  myFab  = mf[mfi];
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();

        for (int kx = fab_lo[2]; kx <= fab_hi[2]; kx++)
        {
            for (int jx = fab_lo[1]; jx <= fab_hi[1]; jx++)
            {
                for (int ix = fab_lo[0]; ix <= fab_hi[0]; ix++)
                {
            	    IntVect indices(D_DECL(ix, jx, kx));

	            for (int n = 0; n < BL_SPACEDIM; n++)
	            {
                        disp[n] = myFab(indices,n);
                        //
			// Start with homogeneous distribution (for 1 p per cell in the center of the cell),
                        // then add the displacement (input values weighted by domain length).
                        //
	                p.m_pos[n] = geom.ProbLo(n) + 
                                      (indices[n]+Real(0.5))*dx[n] +
		                      disp[n] * len[n];
                        //
			// Set the velocities.
                        //
	                p.m_data[n+1] = disp[n] * vel_fac[n];
	            }
                    //
		    // Set the mass of the particle from the input value.
                    //
	            p.m_data[0] = particleMass;
	            p.m_id      = ParticleBase::NextID();
	            p.m_cpu     = MyProc;
	
	            if (!ParticleBase::Where(p,m_gdb))
                    {
      		        ParticleBase::PeriodicShift(p,m_gdb);

                        if (!ParticleBase::Where(p,m_gdb))
                            BoxLib::Abort("ParticleContainer<N>::InitCosmo1ppc(): invalid particle");
		    }

	            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
	            //
	            // Add it to the appropriate PBox at the appropriate level.
	            //
	            m_particles[p.m_lev][p.m_grid].push_back(p);
                }
            }
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass)
{
    Real shift[] = {0,0,0};
    InitCosmo(mf, vel_fac, n_part, particleMass, shift);
}

template <int N>
void
ParticleContainer<N>::InitCosmo  (MultiFab& mf, const Real vel_fac[], const Array<int> n_part, const Real particleMass, const Real shift[])
{
    BL_PROFILE("ParticleContainer<N>::InitCosmo()");
    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    const Real len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };
    //
    // Print the grids as a sanity check.
    //
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const int  *fab_lo = mfi.validbox().loVect();
        const int  *fab_hi = mfi.validbox().hiVect();
        if (mfi.validbox().isEmpty())
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Empty box in InitCosmo ");
        }
        if (!geom.Domain().contains(mfi.validbox()))
        {
           std::cout << "...bad grid lo " << fab_lo[0] << " " << fab_lo[1] << " " << fab_lo[2] << '\n';
           std::cout << "...bad grid hi " << fab_hi[0] << " " << fab_hi[1] << " " << fab_hi[2] << '\n';
           BoxLib::Error("Box in InitCosmo not contained in domain");
        }
    }
    //
    // We will need one ghost cell, so check wether we have one.
    //
    if (mf.nGrow() < 1)
        BoxLib::Abort("ParticleContainer<N>::InitCosmo: mf needs at least one correctly filled ghost zone!");

    if ( !(n_part[0] == n_part[1] && n_part[1] == n_part[2]) )
    {
	    std::cout << '\n' << '\n';
	    std::cout << "Your particle lattice will have different spacings in the spatial directions!" << '\n';
	    std::cout << "You might want to change the particle number or the algorithm... ;)" << '\n';
	    std::cout << '\n' << '\n';
    }
    //
    // Place the particles evenly spaced in the problem domain.
    // Not perfectly fast - but easy
    //
    Real         pos[BL_SPACEDIM];
    ParticleType p;

    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
    {
        const Box&  box     = mfi.validbox();
        RealBox     gridloc = RealBox(box, geom.CellSize(), geom.ProbLo());
	const Real* xlo     = gridloc.lo();
	const Real* xhi     = gridloc.hi();

        for (int k = 0; k < n_part[2]; k++)
        {
	    for (int j = 0; j < n_part[1]; j++)
            {
	        for (int i = 0; i < n_part[0]; i++)
                {
		    bool    isInValidBox = true;
            	    IntVect indices(D_DECL(i, j, k));

		    for (int n = 0; n < BL_SPACEDIM; n++)
                    {
  		        pos[n] = geom.ProbLo(n)
		               + (indices[n] + Real(0.5))*len[n]/n_part[n]
			       + shift[n];
                        //
			// Make sure particle is not on a boundary...
                        //
			pos[n] += 1e-14 * (geom.ProbHi(n) - geom.ProbLo(n));

			isInValidBox = isInValidBox 
				     && (pos[n] > xlo[n]) 
				     && (pos[n] < xhi[n]);
		    }

		    if (isInValidBox)
                    {
                        D_TERM(p.m_pos[0] = pos[0];,
                               p.m_pos[1] = pos[1];,
                               p.m_pos[2] = pos[2];);
                        //
		        // Set the mass of the particle.
                        //
	                p.m_data[0] = particleMass;
	                p.m_id      = ParticleBase::NextID();
	                p.m_cpu     = MyProc;

	                if (!ParticleBase::Where(p,m_gdb))
                        {
      		            ParticleBase::PeriodicShift(p,m_gdb);

                            if (!ParticleBase::Where(p,m_gdb))
                                BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
		        }

	                BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
	                //
	                // Add it to the appropriate PBox at the appropriate level.
	                //
	                m_particles[p.m_lev][p.m_grid].push_back(p);
		    }
	        }
	    }
        }
    }

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with equidistant placement" << '\n';
    }
    //
    // Let Redistribute() sort out where the particles belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Redistribute done" << '\n';
    }

    BL_ASSERT(OK());
    BL_ASSERT(N >= BL_SPACEDIM+1);
    //
    // FIXME: Will we ever need initial particles in grids deeper than 0?!
    //
    PMap& pmap = m_particles[0];
    //
    // Make sure, that mf and m_gdb.boxArray(0) are defined on the same boxarray.
    //
    for (typename PMap::iterator pmap_it = pmap.begin(), pmapEnd = pmap.end(); pmap_it != pmapEnd; ++pmap_it)
    {
        const int        grid    = pmap_it->first;
        PBox&            pbox    = pmap_it->second;
        const int        n       = pbox.size();
        const FArrayBox& dfab    = mf[grid];

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int i = 0; i < n; i++)
        {
            ParticleType& p = pbox[i];

            if (p.m_id <= 0) continue;

            BL_ASSERT(p.m_grid == grid);

            Real disp[BL_SPACEDIM];
            //
	    // Do CIC interpolation onto the particle positions.
	    // For CIC we need one ghost cell!
            //
            ParticleBase::GetGravity(dfab,m_gdb->Geom(p.m_lev),p,disp);

            D_TERM(p.m_pos[0] += len[0]*disp[0];,
                   p.m_pos[1] += len[1]*disp[1];,
                   p.m_pos[2] += len[2]*disp[2];);
            //
            // Note: m_data[0] is mass, 1 is v_x, ...
            //
            D_TERM(p.m_data[1] = vel_fac[0]*disp[0];,
                   p.m_data[2] = vel_fac[1]*disp[1];,
                   p.m_data[3] = vel_fac[2]*disp[2];);

            if (!ParticleBase::Where(p,m_gdb))
            {
	        ParticleBase::PeriodicShift(p,m_gdb);

                if (!ParticleBase::Where(p,m_gdb))
                    BoxLib::Abort("ParticleContainer<N>::InitCosmo(): invalid particle");
	    }

            ParticleBase::Reset(p,m_gdb,true);
        }
    }
    //
    // Let Redistribute() sort out where the particles now belong.
    //
    Redistribute(true);

    if (ParallelDescriptor::IOProcessor() && m_verbose)
    {
        std::cout << "Done with particle displacement" << '\n';
    }

    if (m_verbose > 1)
    {
        Real runtime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(runtime, IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "InitCosmo() done time: " << runtime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitRandom (long          icount,
                                  unsigned long iseed,
                                  Real          mass,
                                  bool          serialize,
                                  RealBox       containing_bx)
{
    BL_PROFILE("ParticleContainer<N>::InitRandom()");
    BL_ASSERT(iseed  > 0);
    BL_ASSERT(icount > 0);

    BL_ASSERT(m_gdb != 0);

    const int       MyProc   = ParallelDescriptor::MyProc();
    const int       NProcs   = ParallelDescriptor::NProcs();
    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    Real r, x, len[BL_SPACEDIM] = { D_DECL(geom.ProbLength(0),
                                           geom.ProbLength(1),
                                           geom.ProbLength(2)) };

    // We will enforce that the particles are within the containing_bx.  
    // If containing_bx is not passed in, it defaults to the full domain.
    if (!containing_bx.ok()) containing_bx = geom.ProbDomain();

    // containing_bx is assumed to lie within the domain.
    BL_ASSERT(geom.ProbDomain().contains(containing_bx));

    const Real* xlo = containing_bx.lo();
    const Real* xhi = containing_bx.hi();

    BoxLib::InitRandom(iseed+MyProc);

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.

    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
    {
        BL_ASSERT(m_particles[lev].empty());
    }

    if (serialize)
    {
        //
        // We'll let IOProc generate the particles so we get the same
        // positions no matter how many CPUs we have.  This is here
        // mainly for debugging purposes.  It's not really useful for
        // very large numbers of particles.
        //
        //
        Array<ParticleBase::RealType> pos(icount*BL_SPACEDIM);

        if (ParallelDescriptor::IOProcessor())
        {
            for (long j = 0; j < icount; j++)
            {
                for (int i = 0; i < BL_SPACEDIM; i++)
                {
                    do
                    {
                        r = BoxLib::Random();
                        x = geom.ProbLo(i) + (r * len[i]);
                    }
                    while (x < xlo[i] || x > xhi[i]);

                    pos[j*BL_SPACEDIM + i] = x;
                }
            }
        }

        ParallelDescriptor::Bcast(pos.dataPtr(), icount*BL_SPACEDIM, IOProc);

        int cnt = 0;

        for (long j = 0; j < icount; j++)
        {
            ParticleType p;

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                p.m_pos[i] = pos[j*BL_SPACEDIM + i];
            }

            if (N > 0)
            {
                p.m_data[0] = mass;

                for (int i = 1; i < N; i++)
                    //
                    // Just zero out the rest of the data for lack of a better value.
                    //
                    p.m_data[i] = 0;
            }

            if (!ParticleBase::Where(p,m_gdb))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());

            const int who = m_gdb->ParticleDistributionMap(p.m_lev)[p.m_grid];

            if (who == MyProc)
            {
                //
                // We own it. Add it to the appropriate PBox at the appropriate level.
                //
                p.m_id  = ParticleBase::NextID();
                p.m_cpu = MyProc;

                m_particles[p.m_lev][p.m_grid].push_back(p);

                cnt++;
            }
        }

        BL_ASSERT(OK());
    }
    else
    {
        //
        // We'll generate the particles in parallel.
        //
        // Each CPU will key off the given seed to get independent streams of random numbers.
        //
        long M = icount / NProcs;
        //
        // Processor 0 will get the slop.
        //
        if (MyProc == 0)
        {
            M += (icount % NProcs);
        }

        for (long icnt = 0; icnt < M; icnt++)
        {
            ParticleType p;

            for (int i = 0; i < BL_SPACEDIM; i++)
            {
                do
                {
                    r = BoxLib::Random();
                    x = geom.ProbLo(i) + (r * len[i]);
                }
                while (x < xlo[i] || x > xhi[i]);

                p.m_pos[i] = x;

                BL_ASSERT(p.m_pos[i] < geom.ProbHi(i));
            }

            if (N > 0)
            {
                p.m_data[0] = mass;

                for (int i = 1; i < N; i++)
                    //
                    // Just zero out the rest of the data for lack of a better value.
                    //
                    p.m_data[i] = 0;
            }

            p.m_id  = ParticleBase::NextID();
            p.m_cpu = MyProc;

            if (!ParticleBase::Where(p,m_gdb))
            {
                BoxLib::Abort("ParticleContainer<N>::InitRandom(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
            //
            // Add it to the appropriate PBox at the appropriate level.
            //
            m_particles[p.m_lev][p.m_grid].push_back(p);
        }
        //
        // Let Redistribute() sort out where the particles belong.
        //
        Redistribute(true);
    }

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::InitRandom() time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitOnePerCell (Real x_off, Real y_off, Real z_off,
                                      Real      mass, 
                                      MultiFab& particle_mf)
{
    BL_PROFILE("ParticleContainer<N>::InitOnePerCell()");

    BL_ASSERT(m_gdb != 0);

    // Note that x_off, y_off and z_off are the offsets from the lower left corner
    //      in each cell as measured in units of dx, so they must be in [0,1]
    BL_ASSERT(x_off >= 0. && y_off >= 0. && z_off >= 0.); 
    BL_ASSERT(x_off <= 1. && y_off <= 1. && z_off <= 1.); 

    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    //
    // This assumes level 0 since geom = m_gdb->Geom(0)
    //
    const Real* dx  = geom.CellSize();

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
        BL_ASSERT(m_particles[lev].empty());

    //
    // We'll generate the particles in parallel -- but no tiling here.
    //
    for (MFIter mfi(particle_mf,false); mfi.isValid(); ++mfi)
    {
        Box grid = m_gdb->ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box = RealBox(grid,dx,geom.ProbLo());

	int nx = grid.length(0), ny = grid.length(1), nz = grid.length(2); 

        for (int i = 0; i < nx; i++) {
         for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++)
          {
            ParticleType p;

            p.m_pos[0] = grid_box.lo(0) + (x_off + i)*dx[0];
            p.m_pos[1] = grid_box.lo(1) + (y_off + j)*dx[1];
            p.m_pos[2] = grid_box.lo(2) + (z_off + k)*dx[2];

            for (int i = 0; i < BL_SPACEDIM; i++)
                BL_ASSERT(p.m_pos[i] < grid_box.hi(i));

            if (N > 0)
            {
                p.m_data[0] = mass;

                for (int i = 1; i < N; i++)
                    //
                    // Just zero out the rest of the data for lack of a better value.
                    //
                    p.m_data[i] = 0;
            }

            p.m_id  = ParticleBase::NextID();
            p.m_cpu = ParallelDescriptor::MyProc();

            if (!ParticleBase::Where(p,m_gdb))
            {
                BoxLib::Abort("ParticleContainer<N>::InitOnePerCell(): invalid particle");
            }

            BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
            //
            // Add it to the appropriate PBox at the appropriate level.
            //
            m_particles[p.m_lev][p.m_grid].push_back(p);
          } 
         } 
        }
    }

    //
    // We still need to redistribute in order to define each particle's cell, grid and level, but this 
    //    shouldn't require any inter-node communication because the particles should already be in the right grid.
    //
    Redistribute(true);

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::InitOnePerCell() time: " << stoptime << '\n';
        }
    }
}

template <int N>
void
ParticleContainer<N>::InitNRandomPerCell (int n_per_cell, Real mass, MultiFab& particle_mf)
{
    BL_PROFILE("ParticleContainer<N>::InitNRandomPerCell()");

    BL_ASSERT(m_gdb != 0);

    const int       IOProc   = ParallelDescriptor::IOProcessorNumber();
    const Real      strttime = ParallelDescriptor::second();
    const Geometry& geom     = m_gdb->Geom(0);

    //
    // This assumes level 0 since geom = m_gdb->Geom(0)
    //
    const Real* dx  = geom.CellSize();

    Real r;

    m_particles.reserve(15);  // So we don't ever have to do any copying on a resize.
    m_particles.resize(m_gdb->finestLevel()+1);

    for (int lev = 0; lev < m_particles.size(); lev++)
        BL_ASSERT(m_particles[lev].empty());

    //
    // We'll generate the particles in parallel -- but no tiling here.
    //
    for (MFIter mfi(particle_mf,false); mfi.isValid(); ++mfi)
    {
        Box grid = m_gdb->ParticleBoxArray(0)[mfi.index()];
        RealBox grid_box = RealBox(grid,dx,geom.ProbLo());

	int nx = grid.length(0), ny = grid.length(1), nz = grid.length(2); 

        for (int i = 0; i < nx; i++) {
         for (int j = 0; j < ny; j++) {
          for (int k = 0; k < nz; k++)
          {
            ParticleType p;

            for (int n = 0; n < n_per_cell; n++)
            {
               r = BoxLib::Random();
               p.m_pos[0] = grid_box.lo(0) + (r + i)*dx[0];
               r = BoxLib::Random();
               p.m_pos[1] = grid_box.lo(1) + (r + j)*dx[1];
               r = BoxLib::Random();
               p.m_pos[2] = grid_box.lo(2) + (r + k)*dx[2];

               for (int i = 0; i < BL_SPACEDIM; i++)
                   BL_ASSERT(p.m_pos[i] < grid_box.hi(i));

               if (N > 0)
               {
                   p.m_data[0] = mass;
   
                   for (int i = 1; i < N; i++)
                       //
                       // Just zero out the rest of the data for lack of a better value.
                       //
                       p.m_data[i] = 0;
               }

               p.m_id  = ParticleBase::NextID();
               p.m_cpu = ParallelDescriptor::MyProc();

               if (!ParticleBase::Where(p,m_gdb))
                   BoxLib::Abort("ParticleContainer<N>::InitNRandomPerCell(): invalid particle");

               BL_ASSERT(p.m_lev >= 0 && p.m_lev <= m_gdb->finestLevel());
               //
               // Add it to the appropriate PBox at the appropriate level.
               //
               m_particles[p.m_lev][p.m_grid].push_back(p);
            }
          } 
         } 
        }
    }

    //
    // We still need to redistribute in order to define each particle's cell, grid and level, but this 
    //    shouldn't require any inter-node communication because the particles should already be in the right grid.
    //
    Redistribute(true);

    if (m_verbose > 1)
    {
        Real stoptime = ParallelDescriptor::second() - strttime;

        ParallelDescriptor::ReduceRealMax(stoptime,IOProc);

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "ParticleContainer<N>::InitNRandomPerCell() time: " << stoptime << '\n';
        }
    }
}
