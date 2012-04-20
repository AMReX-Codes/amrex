#include <Particles.H>
#include <ParmParse.H>
#include <limits>

void
ParticleBase::AssignDensityCoeffs (const ParticleBase& p,
                                   const Real*         plo,
                                   const Real*         dx,
                                   Real*               fracs,
                                   IntVect*            cells)
{
    const Real len[BL_SPACEDIM] = { D_DECL((p.m_pos[0]-plo[0])/dx[0] + 0.5,
                                           (p.m_pos[1]-plo[1])/dx[1] + 0.5,
                                           (p.m_pos[2]-plo[2])/dx[2] + 0.5) };

    const IntVect csect(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    Real frac[BL_SPACEDIM] = { D_DECL(len[0]-csect[0], len[1]-csect[1], len[2]-csect[2]) };

    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (frac[d] > 1) frac[d] = 1;
        if (frac[d] < 0) frac[d] = 0;
    }

    IntVect cell = csect;
    //
    // Note that cells[0] contains csect.
    //
#if (BL_SPACEDIM == 1)
    // High
    cells[0] = cell;
    fracs[0] = frac[0];

    // Low
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]);

#elif (BL_SPACEDIM == 2)
    // HH
    cells[0] = cell;
    fracs[0] = frac[0] * frac[1] ;
    
    // LH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]) * frac[1];
    
    // LL
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    fracs[2] = (1-frac[0]) * (1-frac[1]);
    
    // HL
    cell[0]  = cell[0] + 1;
    cells[3] = cell;
    fracs[3] = frac[0] * (1-frac[1]);

#elif (BL_SPACEDIM == 3)
    // HHH
    cells[0] = cell;
    fracs[0] = frac[0] * frac[1] * frac[2];

    // LHH
    cell[0]  = cell[0] - 1;
    cells[1] = cell;
    fracs[1] = (1-frac[0]) * frac[1] * frac[2];

    // LLH
    cell[1]  = cell[1] - 1;
    cells[2] = cell;
    fracs[2] = (1-frac[0]) * (1-frac[1]) * frac[2];
    
    // HLH
    cell[0]  = cell[0] + 1;
    cells[3] = cell;
    fracs[3] = frac[0] * (1-frac[1]) * frac[2];

    cell = csect;

    // HHL
    cell[2]  = cell[2] - 1;
    cells[4] = cell;
    fracs[4] = frac[0] * frac[1] * (1-frac[2]);
    
    // LHL
    cell[0]  = cell[0] - 1;
    cells[5] = cell;
    fracs[5] = (1-frac[0]) * frac[1] * (1-frac[2]);

    // LLL
    cell[1]  = cell[1] - 1;
    cells[6] = cell;
    fracs[6] = (1-frac[0]) * (1-frac[1]) * (1-frac[2]);
    
    // HLL
    cell[0]  = cell[0] + 1;
    cells[7] = cell;
    fracs[7] = frac[0] * (1-frac[1]) * (1-frac[2]);
#endif
}


bool
ParticleBase::CrseToFine (const BoxArray& cfba,
                          const IntVect*  cells,
                          bool*           which)
{
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle, will we cross a  crse->fine boundary of the level
    // with coarsened fine BoxArray "cfba".  "cells" are as calculated from
    // AssignDensityCoeffs().
    //
    const int M = D_TERM(2,+2,+4);

    for (int i = 0; i < M; i++)
        which[i] =  false;

    bool result = false;

    for (int i = 0; i < M; i++)
    {
        if (cfba.contains(cells[i]))
        {
            result   = true;
            which[i] = true;
        }
    }

    return result;
}

bool
ParticleBase::FineToCrse (const ParticleBase& p,
                          int                 flev,
                          const Amr*          amr,
                          const IntVect*      fcells,
                          const BoxArray&     fvalid,
                          bool*               which,
                          int*                cgrid)
{
    BL_ASSERT(amr != 0);
    BL_ASSERT(flev > 0);
    //
    // We're in AssignDensity(). We want to know whether or not updating
    // with a particle we'll cross a fine->crse boundary.  Note that crossing
    // a periodic boundary, where the periodic shift lies in our valid region,
    // is not considered a Fine->Crse crossing.
    //
    const int M = D_TERM(2,+2,+4);

    for (int i = 0; i < M; i++)
    {
        cgrid[i] = -1;
        which[i] = false;
    }

    const Box ibx = BoxLib::grow(amr->boxArray(flev)[p.m_grid],-1);

    BL_ASSERT(ibx.ok());

    if (ibx.contains(p.m_cell))
        //
        // We're strictly contained in our valid box.
        // We can't cross a fine->crse boundary.
        //
        return false;
    //
    // Otherwise ...
    //
    Real    cfracs[M];
    IntVect ccells[M];

    const Geometry& cgm = amr->Geom(flev-1);

    ParticleBase::AssignDensityCoeffs(p, cgm.ProbLo(), cgm.CellSize(), cfracs, ccells);

    bool result = false;

    std::vector< std::pair<int,Box> > isects;

    Array<IntVect> pshifts(27);

    for (int i = 0; i < M; i++)
    {
        const IntVect iv = ccells[i] * amr->refRatio(flev-1);
        //
        // Note that "fvalid" may not equal the valid BoxArray.  It will
        // also include any ghost regions outside the valid domain, that can be
        // periodically shifted back into the valid region of the BoxArray
        // at this level.
        //
        if (!fvalid.contains(iv))
        {
            result   = true;
            which[i] = true;

            Box cbx(ccells[i],ccells[i]);

            if (!cgm.Domain().contains(ccells[i]))
            {
                //
                // We must be at a periodic boundary.
                // Find valid box into which we can be periodically shifted.
                //
                BL_ASSERT(cgm.isAnyPeriodic());

                cgm.periodicShift(cbx, cgm.Domain(), pshifts);

                BL_ASSERT(pshifts.size() == 1);

                cbx -= pshifts[0];

                BL_ASSERT(cbx.ok());
                BL_ASSERT(cgm.Domain().contains(cbx));
            }
            //
            // Which grid at the crse level do we need to update?
            //
            isects = amr->boxArray(flev-1).intersections(cbx);

            BL_ASSERT(!isects.empty());
            BL_ASSERT(isects.size() == 1);

            cgrid[i] = isects[0].first;  // The grid ID at crse level that we hit.
        }
    }

    return result;
}

void
ParticleBase::FineCellsToUpdateFromCrse (const ParticleBase& p,
                                         int                 lev,
                                         const Amr*          amr,
                                         const IntVect&      ccell,
                                         Array<int>&         fgrid,
                                         Array<Real>&        ffrac,
                                         Array<IntVect>&     fcells)
{
    BL_ASSERT(lev >= 0);
    BL_ASSERT(lev < amr->finestLevel());

    const Box fbx = BoxLib::refine(Box(ccell,ccell),amr->refRatio(lev));

    const BoxArray& fba = amr->boxArray(lev+1);
    const Real*     plo = amr->Geom(lev).ProbLo();
    const Real*     dx  = amr->Geom(lev).CellSize();
    const Real*     fdx = amr->Geom(lev+1).CellSize();

    BL_ASSERT(fba.contains(fbx));

    fgrid.clear();
    ffrac.clear();
    fcells.clear();
    //
    // Which fine cells does particle "p" that wants to update "ccell" do we touch at the finer level?
    //
    for (IntVect iv = fbx.smallEnd(); iv <= fbx.bigEnd(); fbx.next(iv))
    {
        bool touches = true;

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            const Real celllo = iv[k]  * fdx[k] + plo[k];
            const Real cellhi = celllo + fdx[k];

            if ((p.m_pos[k] < celllo) && (celllo > (p.m_pos[k] + dx[k]/2)))
                touches = false;

            if ((p.m_pos[k] > cellhi) && (cellhi < (p.m_pos[k] - dx[k]/2)))
                touches = false;
        }

        if (touches)
            fcells.push_back(iv);
    }

    std::vector< std::pair<int,Box> > isects;

    Real sum_fine = 0;
    //
    // We need to figure out the fine fractions.
    //
    for (int j = 0; j < fcells.size(); j++)
    {
        const IntVect& iv = fcells[j];

        isects = fba.intersections(Box(iv,iv));

        BL_ASSERT(!isects.empty());
        BL_ASSERT(isects.size() == 1);

        fgrid.push_back(isects[0].first);

        Real the_frac = 1;

        for (int k = 0; k < BL_SPACEDIM; k++)
        {
            const Real celllo = (iv[k] * fdx[k] + plo[k]);

            if (p.m_pos[k] <= celllo)
            {
                const Real isecthi = p.m_pos[k] + dx[k]/2;

                the_frac *= std::min((isecthi - celllo),fdx[k]);
            }
            else
            {
                const Real cellhi  = (iv[k]+1) * fdx[k] + plo[k];
                const Real isectlo = p.m_pos[k] - dx[k]/2;

                the_frac *= std::min((cellhi - isectlo),fdx[k]);
            }
        }

        ffrac.push_back(the_frac);

        sum_fine += the_frac;
    }

    BL_ASSERT(ffrac.size() == fcells.size());
    BL_ASSERT(fgrid.size() == fcells.size());
    //
    // Now adjust the fine fractions so they sum to one.
    //
    for (int j = 0; j < ffrac.size(); j++)
    {
        ffrac[j] /= sum_fine;
        if (ffrac[j] > 1)
            ffrac[j] = 1;
    }
}

//
// Used by AssignDensity (PArray<MultiFab>& mf).
//
// Passes data needed by Crse->Fine or Fine->Crse to CPU that needs it.
//
// We store the data that needs to be sent in "data".
//

void
ParticleBase::AssignDensityDoit (PArray<MultiFab>&         mf,
                                 std::deque<ParticleBase>& data)
{
    const int NProcs = ParallelDescriptor::NProcs();

    if (NProcs == 1)
    {
        BL_ASSERT(data.empty());
        return;
    }

#if BL_USE_MPI
    //
    // We may have data that needs to be sent to another CPU.
    //
    const int MyProc = ParallelDescriptor::MyProc();

    Array<int> Snds(NProcs,0);
    Array<int> Rcvs(NProcs,0);

    for (std::deque<ParticleBase>::const_iterator it = data.begin(), End = data.end();
         it != End;
         ++it)
    {
        const int lev = it->m_lev;
        const int grd = it->m_grid;

        BL_ASSERT(lev >= 0 && lev < mf.size());
        BL_ASSERT(grd >= 0 && grd < mf[lev].size());

        const int who = mf[lev].DistributionMap()[grd];

        BL_ASSERT(who != MyProc);
        BL_ASSERT(mf[lev].fabbox(grd).contains(it->m_cell));

        Snds[who]++;
    }

    BL_ASSERT(Snds[MyProc] == 0);

    long maxsendcount = 0;
    for (int i = 0; i < NProcs; i++)
        maxsendcount += Snds[i];
    ParallelDescriptor::ReduceLongMax(maxsendcount);

    if (maxsendcount == 0)
    {
        //
        // There's no parallel work to do.
        //
        BL_ASSERT(data.empty());
        return;
    }

    BL_MPI_REQUIRE( MPI_Alltoall(Snds.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 Rcvs.dataPtr(),
                                 1,
                                 ParallelDescriptor::Mpi_typemap<int>::type(),
                                 ParallelDescriptor::Communicator()) );
    BL_ASSERT(Rcvs[MyProc] == 0);

    int NumRcvs = 0;
    for (int i = 0; i < NProcs; i++)
        NumRcvs += Rcvs[i];

    int NumSnds = 0;
    for (int i = 0; i < NProcs; i++)
        NumSnds += Snds[i];

    BL_ASSERT(data.size() == NumSnds);
    //
    // The data we receive from ParticleBases.
    //
    // We only use: m_lev, m_grid, m_cell & m_pos[0].
    //
    const int iChunkSize = 2 + BL_SPACEDIM;
    const int rChunkSize = 1;

    Array<int>  irecvdata (NumRcvs*iChunkSize);
    Array<Real> rrecvdata (NumRcvs*rChunkSize);

    Array<int>   offset(NProcs);
    Array<int>  sdispls(NProcs);
    Array<int>  rdispls(NProcs);
    Array<int> sendcnts(NProcs);
    Array<int> recvcnts(NProcs);

    {
        //
        // First send/recv "int" data.
        //
        Array<int> senddata (NumSnds*iChunkSize);

        offset[0] = sdispls[0] = rdispls[0] = 0;

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * iChunkSize;
            sendcnts[i] = Snds[i] * iChunkSize;

            if (i > 0)
            {
                offset [i] = offset [i-1] + sendcnts[i-1];
                rdispls[i] = rdispls[i-1] + recvcnts[i-1];
                sdispls[i] = sdispls[i-1] + sendcnts[i-1];
            }
        }

        for (std::deque<ParticleBase>::const_iterator it = data.begin(), End = data.end();
             it != End;
             ++it)
        {
            const int who  = mf[it->m_lev].DistributionMap()[it->m_grid];
            const int ioff = offset[who];

            senddata[ioff+0] = it->m_lev;
            senddata[ioff+1] = it->m_grid;

            D_TERM(senddata[ioff+2] = it->m_cell[0];,
                   senddata[ioff+3] = it->m_cell[1];,
                   senddata[ioff+4] = it->m_cell[2];);

            offset[who] += iChunkSize;
        }

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      NumRcvs == 0 ? 0 : irecvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<int>::type(),
                                      ParallelDescriptor::Communicator()) );
    }

    {
        //
        // Now send/recv the Real data.
        //
        Array<Real> senddata (NumSnds*rChunkSize);

        offset[0] = sdispls[0] = rdispls[0] = 0;

        for (int i = 0; i < NProcs; i++)
        {
            recvcnts[i] = Rcvs[i] * rChunkSize;
            sendcnts[i] = Snds[i] * rChunkSize;

            if (i > 0)
            {
                offset [i] = offset [i-1] + sendcnts[i-1];
                rdispls[i] = rdispls[i-1] + recvcnts[i-1];
                sdispls[i] = sdispls[i-1] + sendcnts[i-1];
            }
        }

        for (std::deque<ParticleBase>::const_iterator it = data.begin(), End = data.end();
             it != End;
             ++it)
        {
            const int who = mf[it->m_lev].DistributionMap()[it->m_grid];

            senddata[offset[who]] = it->m_pos[0];

            offset[who]++;
        }
        //
        // We can free up memory held by "data" -- don't need it anymore.
        //
        std::deque<ParticleBase>().swap(data);

        BL_MPI_REQUIRE( MPI_Alltoallv(NumSnds == 0 ? 0 : senddata.dataPtr(),
                                      sendcnts.dataPtr(),
                                      sdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      NumRcvs == 0 ? 0 : rrecvdata.dataPtr(),
                                      recvcnts.dataPtr(),
                                      rdispls.dataPtr(),
                                      ParallelDescriptor::Mpi_typemap<Real>::type(),
                                      ParallelDescriptor::Communicator()) );
    }
    //
    // Now update "mf".
    //
    if (NumRcvs > 0)
    {
        const int*  idata = irecvdata.dataPtr();
        const Real* rdata = rrecvdata.dataPtr();

        for (int i = 0; i < NumRcvs; i++)
        {
            const int     lev  = idata[0];
            const int     grd  = idata[1];
            const IntVect cell = IntVect(D_DECL(idata[2],idata[3],idata[4]));

            BL_ASSERT(mf[lev].DistributionMap()[grd] == MyProc);
            BL_ASSERT(mf[lev][grd].box().contains(cell));

            mf[lev][grd](cell) += *rdata;

            idata += iChunkSize;

            rdata++;
        }
    }
#endif
}

int
ParticleBase::MaxReaders ()
{
    const int Max_Readers_def = 64;

    static int Max_Readers;

    static bool first = true;

    if (first)
    {
        first = false;

        ParmParse pp("particles");

        Max_Readers = Max_Readers_def;

        pp.query("nreaders", Max_Readers);

        Max_Readers = std::min(ParallelDescriptor::NProcs(),Max_Readers);

        if (Max_Readers <= 0)
            BoxLib::Abort("particles.nreaders must be positive");
    }

    return Max_Readers;
}

const std::string&
ParticleBase::DataPrefix ()
{
    //
    // The actual particle data is stored in files of the form: DATA_nnnn.
    //
    static const std::string data("DATA_");

    return data;
}

const std::string&
ParticleBase::Version ()
{
    //
    // If we change the Checkpoint/Restart format we should increment this.
    //
    static const std::string version("Version_One_Dot_Zero");

    return version;
}

static int the_next_id = 1;

int
ParticleBase::NextID ()
{
    int next;

#ifdef _OPENMP
#pragma omp critical(nextid_lock)
#endif
    {
        if (the_next_id == std::numeric_limits<int>::max())
            BoxLib::Abort("ParticleBase::NextID() -- too many particles");

        next = the_next_id++;
    }

    return next;
}

void
ParticleBase::NextID (int nextid)
{
    the_next_id = nextid;
}

IntVect
ParticleBase::Index (const ParticleBase& p,
                     int                 lev,
                     const Amr*          amr)
{
    BL_ASSERT(amr != 0);
    BL_ASSERT(lev >= 0 && lev <= amr->finestLevel());

    IntVect iv;

    const Geometry& geom = amr->Geom(lev);

    D_TERM(iv[0]=floor((p.m_pos[0]-geom.ProbLo(0))/geom.CellSize(0));,
           iv[1]=floor((p.m_pos[1]-geom.ProbLo(1))/geom.CellSize(1));,
           iv[2]=floor((p.m_pos[2]-geom.ProbLo(2))/geom.CellSize(2)););

    iv += geom.Domain().smallEnd();

    return iv;
}

bool
ParticleBase::Where (ParticleBase& p,
                     const Amr*    amr,
                     bool          update)
{
    BL_ASSERT(amr != 0);

    if (update)
    {
        //
        // We have a valid particle whose position has changed slightly.
        // Try to update m_cell and m_grid smartly.
        //
        BL_ASSERT(p.m_id > 0);
        BL_ASSERT(p.m_grid >= 0 && p.m_grid < amr->boxArray(p.m_lev).size());

        IntVect iv = Index(p,p.m_lev,amr);

        if (p.m_cell == iv)
            //
            // The particle hasn't left its cell.
            //
            return true;

        if (p.m_lev == amr->finestLevel())
        {
            p.m_cell = iv;

            if (amr->boxArray(p.m_lev)[p.m_grid].contains(p.m_cell))
                //
                // It has left its cell but is still in the same grid.
                //
                return true;
        }
    }

    std::vector< std::pair<int,Box> > isects;

    for (int lev = amr->finestLevel(); lev >= 0; lev--)
    {
        IntVect iv = Index(p,lev,amr);

        isects = amr->boxArray(lev).intersections(Box(iv,iv));

        if (!isects.empty())
        {
            BL_ASSERT(isects.size() == 1);

            p.m_lev  = lev;
            p.m_grid = isects[0].first;
            p.m_cell = iv;

            return true;
        }
    }

    return false;
}

void
ParticleBase::PeriodicShift (ParticleBase& p,
                             const Amr*    amr)
{
    //
    // This routine should only be called when Where() returns false.
    //
    BL_ASSERT(amr != 0);
    //
    // We'll use level 0 stuff since ProbLo/ProbHi are the same for every level.
    //
    const Geometry& geom = amr->Geom(0);
    const Box&      dmn  = geom.Domain();
    IntVect         iv   = Index(p,0,amr);
    const Real      eps  = 1.e-13;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        if (!geom.isPeriodic(i)) continue;

        if (iv[i] > dmn.bigEnd(i))
        {
            if (p.m_pos[i] == geom.ProbHi(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_pos[i] += eps;

            p.m_pos[i] -= geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] >= geom.ProbLo(i));
        }
        else if (iv[i] < dmn.smallEnd(i))
        {
            if (p.m_pos[i] == geom.ProbLo(i))
                //
                // Don't let particles lie exactly on the domain face.
                // Force the particle to be outside the domain so the
                // periodic shift will bring it back inside.
                //
                p.m_pos[i] -= eps;

            p.m_pos[i] += geom.ProbLength(i);

            BL_ASSERT(p.m_pos[i] <= geom.ProbHi(i));
        }
    }
    //
    // The particle may still be outside the domain in the case
    // where we aren't periodic on the face out which it travelled.
    //
}

void
ParticleBase::Reset (ParticleBase& p,
                     const Amr*    amr,
                     bool          update)
{
    BL_ASSERT(amr != 0);

    if (!ParticleBase::Where(p,amr,update))
    {
        //
        // Here's where we need to deal with boundary conditions.
        //
        // Attempt to shift the particle back into the domain if it
        // crossed a periodic boundary.  Otherwise (for now) we
        // invalidate the particle.
        //
        ParticleBase::PeriodicShift(p,amr);

        if (!ParticleBase::Where(p,amr))
        {
#ifdef _OPENMP
#pragma omp critical(reset_lock)
#endif
            {
                std::cout << "Invalidating out-of-domain particle: " << p << '\n';
            }

            BL_ASSERT(p.m_id > 0);

            p.m_id = -p.m_id;
        }
    }
}

Real
ParticleBase::InterpDoit (const FArrayBox& fab,
                          const IntVect&   hi,
                          const Real*      frac,
                          int              comp)
{
    IntVect cell = hi;

    Real val = 0;

    const Real ifrac[BL_SPACEDIM] = { D_DECL((1-frac[0]), (1-frac[1]), (1-frac[2])) };

#if (BL_SPACEDIM == 1)
    // High
    val += fab(cell,comp) * frac[0];

    // Low
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0];

#elif (BL_SPACEDIM == 2)
    // HH
    val += fab(cell,comp) *    frac[0]  *    frac[1] ;

    // LH
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1] ;

    // LL
    cell[1]   = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1];

    // HL
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1];

#elif (BL_SPACEDIM == 3)
    // HHH
    val += fab(cell,comp) *    frac[0]  *    frac[1]  *    frac[2] ;
   
    // LHH
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1]  *    frac[2] ;
   
    // LLH
    cell[1] = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1] *    frac[2] ;
   
    // HLH
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1] *    frac[2] ;

    cell    = hi;
    cell[2] = cell[2] - 1;

    // HHL
    val += fab(cell,comp) *    frac[0]  *    frac[1]  *    ifrac[2] ;
   
    // LHL
    cell[0] = cell[0] - 1;
    val += fab(cell,comp) * ifrac[0] *    frac[1]  *    ifrac[2] ;
   
    // LLL
    cell[1] = cell[1] - 1;
    val += fab(cell,comp) * ifrac[0] * ifrac[1] *    ifrac[2] ;
   
    // HLL
    cell[0] = cell[0] + 1;
    val += fab(cell,comp) *    frac[0]  * ifrac[1] *    ifrac[2] ;

#endif

    return val;
}

void
ParticleBase::Interp (const ParticleBase& prt,
                      const Amr*          amr,
                      const FArrayBox&    fab,
                      const int*          idx,
                      Real*               val,
                      int                 cnt)
{
    BL_ASSERT(amr != 0);
    BL_ASSERT(idx != 0);
    BL_ASSERT(val != 0);

    const Geometry& gm = amr->Geom(prt.m_lev);
    const Real*     dx = gm.CellSize();

    const Real len[BL_SPACEDIM] = { D_DECL((prt.m_pos[0]-gm.ProbLo(0))/dx[0] + 0.5,
                                           (prt.m_pos[1]-gm.ProbLo(1))/dx[1] + 0.5,
                                           (prt.m_pos[2]-gm.ProbLo(2))/dx[2] + 0.5) };

    const IntVect hi(D_DECL(floor(len[0]), floor(len[1]), floor(len[2])));

    Real frac[BL_SPACEDIM] = { D_DECL(len[0]-hi[0], len[1]-hi[1], len[2]-hi[2]) };

    for (int d = 0; d < BL_SPACEDIM; d++)
    {
        if (frac[d] > 1) frac[d] = 1;
        if (frac[d] < 0) frac[d] = 0;
    }

    for (int i = 0; i < cnt; i++)
    {
        const int comp = idx[i];

        BL_ASSERT(comp >= 0 && comp < fab.nComp());

        val[i] = ParticleBase::InterpDoit(fab,hi,frac,comp);
    }
}

void
ParticleBase::GetGravity (const FArrayBox&    gfab,
                          const Amr*          amr,
                          const ParticleBase& p,
                          Real*               grav)
{
    BL_ASSERT(amr  != 0);
    BL_ASSERT(grav != 0);

    int idx[BL_SPACEDIM] = { D_DECL(0,1,2) };

    ParticleBase::Interp(p,amr,gfab,idx,grav,BL_SPACEDIM);
}

std::ostream&
operator<< (std::ostream& os, const ParticleBase& p)
{
    os << p.m_id   << ' '
       << p.m_cpu  << ' '
       << p.m_lev  << ' '
       << p.m_grid << ' '
       << p.m_cell << ' ';

    for (int i = 0; i < BL_SPACEDIM; i++)
        os << p.m_pos[i] << ' ';

    if (!os.good())
        BoxLib::Error("operator<<(ostream&,ParticleBase&) failed");

    return os;
}
