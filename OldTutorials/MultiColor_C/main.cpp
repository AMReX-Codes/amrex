
#include <AMReX_MultiFab.H>
#include <AMReX_ParmParse.H>
#include <AMReX_VisMF.H>
#include <AMReX_Geometry.H>
#include <AMReX_BndryData.H>
#include <AMReX_LO_BCTYPES.H>
#include <AMReX_MultiGrid.H>

using namespace amrex;

namespace {
    int ncomp  = 8;
    Real a     = 1.e-3;
    Real b     = 1.0;
    Real sigma = 10.0;  // controls the size of jump
    Real w     = 0.05;  // contols the width of the jump

    int  verbose = 0;
    Real tolerance_rel = 1.e-8;
    Real tolerance_abs = 0.0;
}

extern "C"
{
    void fort_set_rhs(double*, const int*, const int*, int,
		      const double*, double, double, double, double);
    void fort_set_coef(double*, const int*, const int*, 
		       double*, const int*, const int*, 
		       int, const double*, double, double);
}

void setup_rhs(MultiFab& rhs, const Geometry& geom);
void setup_coeffs(MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom);
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp);
void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const Vector<MultiFab*>& beta,
	   const Geometry& geom);
void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const Vector<MultiFab*>& beta, 
		   const Geometry& geom);
void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const Vector<MultiFab*>& beta, 
			    const Geometry& geom);

namespace ParallelContext {
    void init();
}

int main(int argc, char* argv[])
{
    amrex::Initialize(argc,argv);
    ParallelContext::init();

    BoxArray ba;
    Geometry geom;
    {
	ParmParse pp;
	
	pp.query("verbose", verbose);

	int n_cell, max_grid_size;
	pp.get("n_cell", n_cell);
	pp.get("max_grid_size", max_grid_size);

	Box domain(IntVect(AMREX_D_DECL(       0,       0,       0)),
		   IntVect(AMREX_D_DECL(n_cell-1,n_cell-1,n_cell-1)));

	ba.define(domain);
	ba.maxSize(max_grid_size);

	RealBox real_box;
	for (int n = 0; n < BL_SPACEDIM; n++) {
	    real_box.setLo(n, 0.0);
	    real_box.setHi(n, 1.0);
	}

	int coord = 0;
 
	geom.define(domain,&real_box,coord);
    }

    DistributionMapping dm{ba};

    MultiFab rhs(ba, dm, ncomp, 0);
    setup_rhs(rhs, geom);

    MultiFab alpha(ba, dm, ncomp, 0);
    Vector<std::unique_ptr<MultiFab> > beta(BL_SPACEDIM);
    for (int i = 0; i < BL_SPACEDIM; ++i) {
        beta[i].reset(new MultiFab(amrex::convert(ba, IntVect::TheDimensionVector(i)),
                                   dm, ncomp, 0));
    }
    setup_coeffs(alpha, amrex::GetVecOfPtrs(beta), geom);

    MultiFab soln(ba, dm, ncomp, 0);

    solve(soln, rhs, alpha, amrex::GetVecOfPtrs(beta), geom);

    VisMF::Write(soln, "soln");

    amrex::Finalize();
}

void setup_rhs(MultiFab& rhs, const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    for ( MFIter mfi(rhs); mfi.isValid(); ++mfi )
    {
	const int* rlo = rhs[mfi].loVect();
	const int* rhi = rhs[mfi].hiVect();

	fort_set_rhs(rhs[mfi].dataPtr(),rlo, rhi, rhs.nComp(),
		     dx, a, b, sigma, w);
    }
}

void setup_coeffs(MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom)
{
    const Real* dx = geom.CellSize();

    alpha.setVal(1.0);

#if (BL_SPACEDIM == 3)
    amrex::Abort("2D only");
#endif

    for ( MFIter mfi(alpha); mfi.isValid(); ++mfi ) 
    {
	FArrayBox& betax = (*beta[0])[mfi];
	FArrayBox& betay = (*beta[1])[mfi];

	fort_set_coef(betax.dataPtr(), betax.loVect(), betax.hiVect(),
		      betay.dataPtr(), betay.loVect(), betay.hiVect(),
		      beta[0]->nComp(), dx, sigma, w);
    }
}

// Dirichlet only in this test
void set_boundary(BndryData& bd, const MultiFab& rhs, int comp)
{
    Real bc_value = 0.0;

    const Geometry& geom = bd.getGeom();
    const Real* dx = geom.CellSize();

    for (int n=0; n<BL_SPACEDIM; ++n) {
	for (MFIter mfi(rhs); mfi.isValid(); ++mfi ) {
	    int i = mfi.index(); 
	    
	    const Box& bx = mfi.validbox();
	    
	    // Our default will be that the face of this grid is either touching another grid
	    //  across an interior boundary or a periodic boundary.  We will test for the other
	    //  cases below.
	    {
		// Define the type of boundary conditions to be Dirichlet (even for periodic)
		bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
		bd.setBoundCond(Orientation(n, Orientation::high),i,comp,LO_DIRICHLET);
	
		// Set the boundary conditions to the cell centers outside the domain
		bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.5*dx[n]);
		bd.setBoundLoc(Orientation(n, Orientation::high),i,0.5*dx[n]);
	    }
	    
	    // Now test to see if we should override the above with Dirichlet physical bc's

	    // We are on the low side of the domain in coordinate direction n
	    if (bx.smallEnd(n) == geom.Domain().smallEnd(n)) {
		// Set the boundary conditions to live exactly on the faces of the domain
		bd.setBoundLoc(Orientation(n, Orientation::low) ,i,0.0 );
	  
		// Set the Dirichlet/Neumann boundary values 
		bd.setValue(Orientation(n, Orientation::low) ,i, bc_value);
	  
		// Define the type of boundary conditions 
		bd.setBoundCond(Orientation(n, Orientation::low) ,i,comp,LO_DIRICHLET);
	    }
	
	    // We are on the high side of the domain in coordinate direction n
	    if (bx.bigEnd(n) == geom.Domain().bigEnd(n)) {
		// Set the boundary conditions to live exactly on the faces of the domain
		bd.setBoundLoc(Orientation(n, Orientation::high) ,i,0.0 );
		
		// Set the Dirichlet/Neumann boundary values
		bd.setValue(Orientation(n, Orientation::high) ,i, bc_value);
		
		// Define the type of boundary conditions 
		bd.setBoundCond(Orientation(n, Orientation::high) ,i,comp,LO_DIRICHLET);
	    }
	}
    }
}

namespace {
    template <class T>
    T sum(Vector<T> v) {
        T s = 0;
        for (int i = 0; i < v.size(); ++i) {
            s += v[i];
        }
        return s;
    }
}

namespace ParallelContext {

    struct Frame
    {
        MPI_Comm comm;

#if 0
        // potential implementation where members of task may have non-contiguous global ranks

        // ranks that are members of this frame
        // indexed by rank in current frame
        Vector<int> global_ranks;
        int rank_n() { return global_ranks.size(); }
        int global_rank(int rank) {
            assert(rank < ranks.size();
            return ranks[rank];
        }
#else
        // members of task must have contiguous global ranks
        int glo_rank_lo, glo_rank_hi;
        int rank_n() { return glo_rank_hi - glo_rank_lo; }
        int glo_rank(int rank) {
            assert(rank >= 0 && rank < rank_n());
            return glo_rank_lo + rank;
        }
#endif

        // specific to this rank
        int loc_rank_me; // local rank relative to current sub-communicator
        int rank_me() { return loc_rank_me; }

        Frame(const Frame &) = delete;
        Frame(Frame &&) = default;
        Frame(MPI_Comm c, int glo, int ghi, int lme) :
          comm(c), glo_rank_lo(glo), glo_rank_hi(ghi), loc_rank_me(lme) { }
        ~Frame() = default;
        Frame &operator=(const Frame &) = delete;
        Frame &operator=(Frame &&) = delete;
    };

    Vector<Frame> frames; // stack of communicator frames

    // call at beginning of program, after MPI_Init
    // probably somewhere inside amrex::Initialize()
    void init() {
        // initialize first frame in stack to global communicator
        int glo_rank_n, glo_rank_me;
        MPI_Comm_size(MPI_COMM_WORLD, &glo_rank_n);
        MPI_Comm_rank(MPI_COMM_WORLD, &glo_rank_me);
        frames.emplace_back(MPI_COMM_WORLD, 0, glo_rank_n, glo_rank_me);
    }

    // number of ranks in world communicator
    int rank_n_global() { return frames[0].rank_n(); }
    // my rank in world communicator
    int rank_me_global() { return frames[0].rank_me(); }
    MPI_Comm comm_global() { return frames[0].comm; }

    // number of ranks in current frame
    int rank_n() { return frames.back().rank_n(); }
    // my sub-rank in current frame
    int rank_me() { return frames.back().rank_me(); }
    MPI_Comm comm() { return frames.back().comm; }

    // split top frame of stack and push new frame on top
    int split(Vector<int> task_rank_n)
    {
        auto task_n = task_rank_n.size();
        assert(sum(task_rank_n) == rank_n());

        // figure out what color (task_me) to pass into MPI_Comm_split
        int glo_rank_lo, glo_rank_hi, loc_rank_me;
        int task_me;
        int lo = 0, hi = 0; // local rank cursors
        for (task_me = 0; task_me < task_n; ++task_me) {
            hi += task_rank_n[task_me];
            if (rank_me() >= lo && rank_me() < hi) {
                glo_rank_lo = glo_rank_lo + lo;
                glo_rank_hi = glo_rank_lo + hi;
                loc_rank_me = rank_me() - lo;
                break;
            }
            lo = hi;
        }
        assert(task_me < task_n);

        MPI_Comm new_comm;
        MPI_Comm_split(comm(), task_me, rank_me(), &new_comm);

        frames.emplace_back(new_comm, glo_rank_lo, glo_rank_hi, loc_rank_me);
        return task_me;
    }

    void unsplit() {
        MPI_Comm_free(&frames.back().comm);
        frames.pop_back();
    }
}

class ForkJoin
{
  public:

    enum Strategy {
        DUPLICATE,  // all tasks get a copy of whole MF
        SPLIT,      // split MF components across tasks
    };
    enum Access { RW, RD, WR };

    struct MFFork
    {
        MultiFab *orig;
        Strategy strategy;
        Access access;
        Vector<MultiFab> forked; // holds new multifab for each task in fork

        MFFork() = default;
        MFFork(MultiFab *o, Strategy s = DUPLICATE, Access a = RW) : orig(o), strategy(s), access(a) {}
    };

    ForkJoin(Vector<int> trn) : task_rank_n(trn) {
        auto rank_n = ParallelContext::rank_n(); // number of ranks in current frame
        auto task_n = task_rank_n.size();
        assert(task_n >= 2);
        assert(sum(task_rank_n) == rank_n);
    }

    ForkJoin(Vector<double> task_rank_pct) {
        auto rank_n = ParallelContext::rank_n(); // number of ranks in current frame
        auto task_n = task_rank_pct.size();
        assert(task_n >= 2);
        task_rank_n.resize(task_n);
        int prev = 0;
        double accum = 0;
        for (int i = 0; i < task_n; ++i) {
            accum += task_rank_pct[i];
            int cur = std::round(task_n * accum);
            task_rank_n[i] = cur - prev;
            prev = cur;
        }
        assert(sum(task_rank_n) == rank_n);
    }

    ForkJoin(int task_n = 2) {
        auto rank_n = ParallelContext::rank_n(); // number of ranks in current frame
        assert(task_n >= 2);
        task_rank_n.resize(task_n);
        for (int i = 0; i < task_n; ++i) {
            task_rank_n[i] = rank_n * (i+1) / task_n -
                             rank_n *  i    / task_n;
        }
    }

    void reg_mf(MultiFab &mf, std::string name, int idx, Strategy strategy, Access access) {
        if (idx >= data[name].size()) {
            data[name].resize(idx + 1);
        }
        data[name][idx] = MFFork(&mf, strategy, access);
    };

    void reg_mf(MultiFab &mf, std::string name, Strategy strategy, Access access) {
        reg_mf(mf, name, 0, strategy, access);
    };

    void reg_mf_vec(Vector<MultiFab *> &mfs, std::string name, Strategy strategy, Access access) {
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, access);
        }
    }

    // these overloads are for in case the MultiFab argument is const
    // access must be read-only
    void reg_mf(const MultiFab &mf, std::string name, int idx, Strategy strategy, Access access) {
        assert(access == RD);
        if (idx >= data[name].size()) {
            data[name].resize(idx + 1);
        }
        data[name][idx] = MFFork(const_cast<MultiFab *>(&mf), strategy, access);
    };

    void reg_mf(const MultiFab &mf, std::string name, Strategy strategy, Access access) {
        reg_mf(mf, name, 0, strategy, access);
    };

    void reg_mf_vec(const Vector<MultiFab *> &mfs, std::string name, Strategy strategy, Access access) {
        for (int i = 0; i < mfs.size(); ++i) {
            reg_mf(*mfs[i], name, i, strategy, access);
        }
    }

    void copy_data_to_tasks() {
        for (auto &p : data) { // for each name
            for (auto &mff : p.second) { // for each index
                const MultiFab &orig = *mff.orig;
                Vector<MultiFab> &forked = mff.forked;
                assert(forked.size() == 0); // should be initialized to empty

                int task_n = task_rank_n.size();
                int comp_n = orig.nComp(); // number of components in original

                forked.reserve(task_n);
                for (int i = 0; i < task_n; ++i) {

                    // create distribution map of current box array over current task's ranks
                    // TODO: hard coded colors only right now
                    //       replace with distribution mapping over custom set of ranks

                    assert(task_n == ParallelDescriptor::NColors());
                    ParallelDescriptor::Color color = ParallelDescriptor::Color(i);
                    int nprocs = ParallelDescriptor::NProcs();
                    const auto &ba = orig.boxArray();
                    DistributionMapping dm {ba, nprocs, color};

                    int comp_lo, comp_hi;
                    if (mff.strategy == SPLIT) {
                        // split components across tasks
                        comp_lo = comp_n *  i    / task_n;
                        comp_hi = comp_n * (i+1) / task_n;
                    } else {
                        // duplicate all components to all tasks
                        comp_lo = 0;
                        comp_hi = comp_n;
                    }

                    // build current task's multifab
                    MultiFab mf(ba, dm, comp_hi - comp_lo, 0);

                    if (mff.access == RD || mff.access == RW) {
                        // parallel copy data into current task's multifab
                        mf.copy(orig, comp_lo, 0, comp_hi - comp_lo);
                    }

                    forked.push_back(std::move(mf));
                }
            }
        }
    }

    void copy_data_from_tasks() {
        for (auto &p : data) { // for each name
            for (auto &mff : p.second) { // for each index
                if (mff.access == WR || mff.access == RW) {
                    MultiFab &orig = *mff.orig;
                    int comp_n = orig.nComp(); // number of components in original
                    const Vector<MultiFab> &forked = mff.forked;
                    if (mff.strategy == SPLIT) {
                        // gather components from across tasks
                        auto task_n = task_rank_n.size();
                        for (int i = 0; i < task_n; ++i) {
                            int comp_lo = comp_n *  i    / task_n;
                            int comp_hi = comp_n * (i+1) / task_n;
                            orig.copy(forked[i], 0, comp_lo, comp_hi - comp_lo);
                        }
                    } else { // mff.strategy == DUPLICATE
                        // TODO: maybe user wants to specify which task's copy to keep?
                        // for now, copy all components from task 0
                        int output_task = 0;
                        orig.copy(forked[output_task], 0, 0, comp_n);
                    }
                }
            }
        }
    }

    template <class F>
    void fork_join(const F &fn)
    {
        copy_data_to_tasks(); // move data to local tasks
        task_me = ParallelContext::split(task_rank_n);
        fn(*this);
        ParallelContext::unsplit();
        copy_data_from_tasks(); // move local data back
    }

    MultiFab &get_mf(std::string name, int idx = 0) {
        assert(idx < data.at(name).size() &&
               task_me < data.at(name)[idx].forked.size());
        return data[name][idx].forked[task_me];
    }

    // vector of pointers to all MFs under a name
    Vector<MultiFab *> get_mf_vec(std::string name) {
        int dim = data.at(name).size();
        Vector<MultiFab *> result(dim);
        for (int idx = 0; idx < dim; ++idx) {
            result[idx] = &data[name][idx].forked[task_me];
        }
        return result;
    }

  private:
    Vector<int> task_rank_n; // number of ranks in each forked task
    int task_me; // task the current rank belongs to
    std::unordered_map<std::string, Vector<MFFork>> data;
};

void solve(MultiFab& soln, const MultiFab& rhs, 
	   const MultiFab& alpha, const Vector<MultiFab*>& beta, const Geometry& geom)
{
    ForkJoin fj(ParallelDescriptor::NColors()); // even split ranks among NColors() tasks
    fj.reg_mf(rhs, "rhs", ForkJoin::SPLIT, ForkJoin::RD);
    fj.reg_mf(alpha, "alpha", ForkJoin::SPLIT, ForkJoin::RD);
    fj.reg_mf_vec(beta, "beta", ForkJoin::SPLIT, ForkJoin::RD);
    // soln gets initialized to 0.0 in single_component_solve
    fj.reg_mf(soln, "soln", ForkJoin::SPLIT, ForkJoin::WR);

    fj.fork_join([&] (ForkJoin &f) {
        colored_solve(f.get_mf("soln"),
                      f.get_mf("rhs"),
                      f.get_mf("alpha"),
                      f.get_mf_vec("beta"),
                      geom);
    });
}

void colored_solve(MultiFab& soln, const MultiFab& rhs, 
		   const MultiFab& alpha, const Vector<MultiFab*>& beta, 
		   const Geometry& geom)
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();

    if (rhs.nComp() == 1)
    {
	single_component_solve(soln, rhs, alpha, beta, geom);
    }
    else
    {	    
	for (int i = 0; i < soln.nComp(); ++i) {
	    MultiFab ssoln (ba, dm, 1, 1);
	    MultiFab srhs  (ba, dm, 1, 0);
	    MultiFab salpha(ba, dm, 1, 0);
	    Vector<std::unique_ptr<MultiFab> > sbeta(BL_SPACEDIM);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		sbeta[j].reset(new MultiFab(beta[j]->boxArray(), dm, 1, 0));
	    }
	    
	    MultiFab::Copy(ssoln , soln , i, 0, 1, 0);
	    MultiFab::Copy(srhs  , rhs  , i, 0, 1, 0);
	    MultiFab::Copy(salpha, alpha, i, 0, 1, 0);
	    for (int j = 0; j < BL_SPACEDIM; ++j) {
		MultiFab::Copy(*sbeta[j], *beta[j], i, 0, 1, 0);
	    }
	    
	    single_component_solve(ssoln, srhs, salpha,
				   amrex::GetVecOfPtrs(sbeta), geom);
	    
	    MultiFab::Copy(soln, ssoln, 0, i, 1, 0);
	}
    }
}

void single_component_solve(MultiFab& soln, const MultiFab& rhs, 
			    const MultiFab& alpha, const Vector<MultiFab*>& beta, 
			    const Geometry& geom)
{
    const BoxArray& ba = soln.boxArray();
    const DistributionMapping& dm = soln.DistributionMap();
    const Real* dx = geom.CellSize();

    BndryData bd(ba, dm, 1, geom);
    set_boundary(bd, rhs, 0);

    ABecLaplacian abec_operator(bd, dx);
    abec_operator.setScalars(a, b);
    abec_operator.setCoefficients(alpha, beta);

    MultiGrid mg(abec_operator);
    mg.setVerbose(verbose);

    soln.setVal(0.0);

    mg.solve(soln, rhs, tolerance_rel, tolerance_abs);
}

