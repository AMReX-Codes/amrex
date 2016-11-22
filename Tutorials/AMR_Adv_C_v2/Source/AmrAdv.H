#ifndef AMR_Adv_C_V2_Adv_H_
#define AMR_Adv_C_V2_Adv_H_

#include <string>
#include <limits>
#include <memory>

#ifdef _OPENMP
#include <omp.h>
#endif

#include <AmrCore.H>
#include <FluxRegister.H>

class AmrAdv
    : public AmrCore
{
public:
    AmrAdv ();
    virtual ~AmrAdv();

    void InitData ();
    void Evolve ();

protected:
    virtual void ErrorEst (int lev, TagBoxArray& tags, Real time, int /*ngrow*/) override;

private:

    void ReadParameters();

    void MakeNewLevel (int lev, Real time,
		       const BoxArray& new_grids, const DistributionMapping& new_dmap);
    void RemakeLevel (int lev, Real time,
		      const BoxArray& new_grids, const DistributionMapping& new_dmap);
    void ClearLevel (int lev);
    void regrid (int lbase, Real time);

    void InitFromScratch ();
    void InitFromCheckpoint ();
    void InitLevelData (int lev);

    void AverageDown ();
    void AverageDownTo (int crse_lev);
    long CountCells (int lev);

    void ComputeDt ();
    Real EstTimeStep (int lev, bool local=false) const;
    void timeStep (int lev, Real time, int iteration);
    void Advance (int lev, Real time, Real dt, int iteration, int ncycle);

    void FillPatch (int lev, Real time, MultiFab& mf, int icomp, int ncomp);
    void GetData (int lev, Real time, PArray<MultiFab>& data, std::vector<Real>& datatime);

    std::string PlotFileName (int lev) const;
    Array<const MultiFab*> PlotFileMF () const;
    Array<std::string> PlotFileVarNames () const;
    void WritePlotFile () const;

    Array<int> istep;      // which step?
    Array<int> nsubsteps;  // how many substeps on each level?
    Array<int> last_regrid_step;  // last time the grid was built

    Array<Real> t_new;
    Array<Real> t_old;
    Array<Real> dt;

    Array<std::unique_ptr<MultiFab> > phi_new;
    Array<std::unique_ptr<MultiFab> > phi_old;

    Array<std::unique_ptr<FluxRegister> > flux_reg; 
    
    // runtime parameters

    int max_step   = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    Real cfl = 0.7;

    int regrid_int = 2;

    int do_reflux = 1;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;
};

#endif
