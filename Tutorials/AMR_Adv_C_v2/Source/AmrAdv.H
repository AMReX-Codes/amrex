#ifndef AMR_Adv_C_V2_Adv_H_
#define AMR_Adv_C_V2_Adv_H_

#include <string>
#include <limits>
#include <memory>

#include <AmrCore.H>

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

    void InitFromScratch ();
    void InitFromCheckpoint ();
    void InitLevelData (int lev);

    void AverageDown ();

    std::string PlotFileName (int lev) const;
    Array<const MultiFab*> PlotFileMF () const;
    std::vector<std::string> PlotFileVarNames () const;

    Array<int> istep;      // which step?
    Array<int> isubstep;   // which substep in subcycle?
    Array<int> nsubsteps;  // how many substeps on each level?

    Array<Real> t_new;
    Array<Real> t_old;
    Array<Real> dt;
    Array<Real> dt_min;

    Array<std::unique_ptr<MultiFab> > phi_new;
    Array<std::unique_ptr<MultiFab> > phi_old;
    
    // runtime parameters

    int max_step   = std::numeric_limits<int>::max();
    Real stop_time = std::numeric_limits<Real>::max();

    Real cfl = 0.7;

    int regrid_int = 2;

    std::string restart_chkfile;

    std::string check_file {"chk"};
    std::string plot_file {"plt"};
    int check_int = -1;
    int plot_int = -1;
};

#endif
