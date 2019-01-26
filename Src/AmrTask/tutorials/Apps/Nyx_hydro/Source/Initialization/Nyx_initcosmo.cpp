#include <iomanip>

#ifdef GRAVITY
#include <Gravity.H>
extern "C"
{ void fort_get_grav_const(amrex::Real* Gconst); }
#endif

#include <Nyx.H>
#include <Nyx_F.H>
#include <NyxParticleContainer.H>
#include <AMReX_Particles_F.H>

using namespace amrex;

#ifdef GRAVITY
void Nyx::icReadAndPrepareFab(std::string mfDirName, int nghost, MultiFab &mf)
{
    if (level > 0 && nghost > 0)
    {
       std::cout << "Are sure you want to do what you are doing?" << std::endl;
       amrex::Abort();
    }

    //
    // Read from the directory into mf
    //
    MultiFab mf_read;

    if (!mfDirName.empty() && mfDirName[mfDirName.length()-1] != '/')
       mfDirName += '/';
    std::string Level = amrex::Concatenate("Level_", level, 1);
    mfDirName.append(Level);
    mfDirName.append("/Cell");

    VisMF::Read(mf_read,mfDirName.c_str());

    if (ParallelDescriptor::IOProcessor())
        std::cout << "mf read" << '\n';

    if (mf_read.contains_nan())
    {
        for (int i = 0; i < mf_read.nComp(); i++)
        {
            if (mf_read.contains_nan(i, 1))
            {
                std::cout << "Found NaNs in read_mf in component " << i << ". " << std::endl;
                amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
            }
        }
    }

    const auto& ba      = parent->boxArray(level);
    const auto& dm      = parent->DistributionMap(level);
    const auto& ba_read = mf_read.boxArray();
    int      nc      = mf_read.nComp();

    //if we don't use a cic scheme for the initial conditions, 
    //we can safely set the number of ghost cells to 0
    //for multilevel ICs we can't use ghostcells
    mf.define(ba, dm, nc, nghost);

    mf.copy(mf_read,0,0,nc);

    mf_read.clear();

    if (! ((ba.contains(ba_read) && ba_read.contains(ba))) )
    {
	if (ParallelDescriptor::IOProcessor()){
            std::cout << "ba      :" << ba << std::endl;
	    std::cout << "ba_read :" << ba_read << std::endl;
            std::cout << "Read mf and hydro mf DO NOT cover the same domain!" 
        	      << std::endl;
	}
	ParallelDescriptor::Barrier();
	if (ParallelDescriptor::IOProcessor()){
            amrex::Abort();
	}
    }

    mf.FillBoundary();
    mf.EnforcePeriodicity(geom.periodicity());

    //FIXME
    //mf.setVal(0);

    if (mf.contains_nan())
    {
        for (int i = 0; i < mf.nComp(); i++)
        {
            if (mf.contains_nan(i, 1, nghost))
            {
                std::cout << "Found NaNs in component " << i << ". " << std::endl;
                amrex::Abort("Nyx::init_particles: Your initial conditions contain NaNs!");
            }
        }
    }
}
#endif

void Nyx::initcosmo()
{

//     if(parent->useFixedUpToLevel()<level)
//       return;
    if (ParallelDescriptor::IOProcessor())
       std::cout << "Calling InitCosmo for level " << level << std::endl;

#ifdef GRAVITY
    Real comoving_OmL;
    Real Gconst;
    const Real len[BL_SPACEDIM] = {geom.ProbLength(0),geom.ProbLength(1),geom.ProbLength(2)};

    Real particleMass;
    std::string mfDirName;

    Real redshift=-1;
    Vector<int> n_part(BL_SPACEDIM);

    if (level > parent->useFixedUpToLevel())
    {
        std::cout << "You have more refinement than grids, there might be a problem with your refinement criterion..." << std::endl;

    	MultiFab& S_new = get_level(level).get_new_data(State_Type);
    	MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

        FillCoarsePatch(S_new, 0, 0,   State_Type, 0, S_new.nComp());
        FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
	return;
    }

    //
    // Read the init directory name and particle mass from the inputs file
    //
    ParmParse pp("cosmo");
    pp.get("initDirName", mfDirName);
#ifdef NUFLUID
    std::string nuMfDirName;
    pp.get("nuInitDirName", nuMfDirName);
#endif
    ParmParse pp2("nyx");
    pp2.get("initial_z",redshift);
    pp2.getarr("n_particles",n_part,0,BL_SPACEDIM);
    
#ifdef NUFLUID
    Real comoving_OmNu;
    fort_get_omnu(&comoving_OmNu);
#endif
    // fort_get_omm(&comoving_OmM );
    // fort_get_omb(&comoving_OmB );
    // fort_get_hubble(&comoving_h);
    fort_get_grav_const(&Gconst);

    // We now define this here instead of reading it.
    comoving_OmL = 1. - comoving_OmM;

    // //For whatever reason OmB is divided by OmM, which I need to undo here...
    // Real realOmB = comoving_OmB*comoving_OmM;

    //compute \Omega_K - just for checking flatness...
    Real comoving_OmK = 1 - comoving_OmM - comoving_OmL;
    if ( std::abs(comoving_OmK) > 1e-3 )
    {
            if (ParallelDescriptor::IOProcessor())
            {
        	    std::cout << "You assume a non-flat universe - \\Omega_K = "
        		      << comoving_OmK << std::endl;
            }
            amrex::Abort();
    }

    //compute \rho_{baryon}=3H_0^2\Omega_{baryon}/(8\pi G)
    Real rhoB = 3*comoving_h*100*comoving_h*100*comoving_OmB 
                / (8*M_PI*Gconst);
    if (ParallelDescriptor::IOProcessor())
    {
       std::cout << "Mean baryonic matter density is " << rhoB 
      	         << " M_sun/Mpc^3." << '\n';
    }

#ifdef NUFLUID
    //compute \rho_{\nu}=3H_0^2\Omega_{\nu}/(8\pi G)
    Real rhoNu = 3*comoving_h*100*comoving_h*100*comoving_OmNu
                 / (8*M_PI*Gconst);
    if(ParallelDescriptor::IOProcessor()){
       std::cout << "Mean neutrino matter density is " << rhoNu
                 << " M_sun/Mpc^3." << std::endl;
    }
#endif

    //compute \rho_{dark matter}=3H_0^2\Omega_{dark matter}/(8\pi G)
#ifdef NUFLUID
    Real comoving_OmD = comoving_OmM - comoving_OmB - comoving_OmNu;
#else
    Real comoving_OmD = comoving_OmM - comoving_OmB;
#endif
    Real rhoD = 3*comoving_h*100*comoving_h*100*comoving_OmD 
                / (8*M_PI*Gconst);
    if (ParallelDescriptor::IOProcessor())
    {
       std::cout << "Mean dark matter density is " << rhoD 
                  << " M_sun/Mpc^3." << '\n';
    }

    if (!do_hydro && comoving_OmB>0){
       if (ParallelDescriptor::IOProcessor()){
          std::cout << std::endl;
          std::cout << std::endl;
          std::cout << "You chose a baryonic matter content greater than 0, "
        	    << "but chose not to do hydro" << std::endl;
          std::cout << "I will be using \\rho_M for the dark matter density..." << std::endl;
          std::cout << std::endl;
          std::cout << std::endl;
       }    
       rhoD += rhoB;
    } 

    //Reads the mf and checks the data...
    MultiFab mf;
    icReadAndPrepareFab(mfDirName, 0, mf);
#ifdef NUFLUID
    MultiFab nuMf;
    icReadAndPrepareFab(nuMfDirName, 0, nuMf);
#endif

    // we have to calculate the initial a on our own
    // as there is no code path to get_comoving_a 
    // (Castro.H sources Particles.H)
    Real comoving_a = 1/(1+redshift);

    
    std::string icSource;
    pp.get("ic-source", icSource);
    int baryon_den, baryon_vx;
    int part_dx, part_vx;
    Real vel_fac[BL_SPACEDIM], dis_fac[BL_SPACEDIM];
    const Real* dx = geom.CellSize();

#ifdef NUFLUID
    int nu_den, nu_vx, nu_vy, nu_vz;
#endif

    if (icSource == "MUSIC")
    {
       if (do_hydro)
       {
          baryon_den = 0;
          baryon_vx = 1;
          part_dx = 4;
          part_vx = 7;
#ifdef NUFLUID
          nu_den = -1;
          nu_vx = -1;
#endif
       }
       else
       {
          baryon_den = -1;
          baryon_vx = -1;
          part_dx = 0;
          part_vx = 3;
#ifdef NUFLUID
          nu_den = -1;
          nu_vx = -1;
#endif
       }
       for (int n=0; n < BL_SPACEDIM; n++)
       {
          //vel_fac[n] = comoving_a * len[n]/comoving_h;
          //vel_fac[n] = len[n]/comoving_h;
          vel_fac[n] = len[n] * comoving_h; // we need the boxlength in Mpc/h
          dis_fac[n] = len[n];
       }
       particleMass = rhoD * dx[0] * dx[1] * dx[2];
       if (ParallelDescriptor::IOProcessor())
       {
          std::cout << "Particle mass at level " << level 
		    << " is " << particleMass 
                    << " M_sun." << '\n';
       }
    }
    else if (icSource == "nyx")
    {
       baryon_den = 3;
       baryon_vx = 0;
       part_dx = 0;
       part_vx = 0;
#ifdef NUFLUID
       nu_den = 3;
       nu_vx = 0;
#endif
       // Velocities are proportional to displacements by
       // LBox [MPc] * a * H [km/s/MPc]
       for (int n=0;n<BL_SPACEDIM;n++)
       {
          vel_fac[n] = len[n]*comoving_a*std::sqrt(comoving_OmM/pow(comoving_a,3)+comoving_OmL)*comoving_h*100;
	  dis_fac[n] = len[n];
       }
       //Compute particle mass
       Real simulationVolume  = len[0]*len[1]*len[2];
       Real  numberOfParticles = n_part[0] * n_part[1] * n_part[2];
       particleMass = rhoD * simulationVolume / numberOfParticles;
       if (ParallelDescriptor::IOProcessor())
       {
          std::cout << "Particle mass is " << particleMass 
                    << " M_sun." << '\n';
       }
    }
    else
    {
       std::cout << "No clue from which code the initial coniditions originate..." << std::endl
	         << "Aborting...!" << std::endl;
       amrex::Abort();
    }

    
    BoxArray baWhereNot;
    if (level < parent->initialBaLevels())
    {
//       baWhereNot = parent->boxArray(level+1);
//       baWhereNot.coarsen(parent->refRatio()[level]);
        baWhereNot = parent->initialBa(level+1);
	//std::cout << "Don't generate particles in region: " << baWhereNot << std::endl;
    }
    //initialize particles
    //std::cout << "FIXME!!!!!!!!!!!!!!!!!!" << std::endl;
    BoxArray myBaWhereNot(baWhereNot.size());
    for (int i=0; i < baWhereNot.size(); i++)
       myBaWhereNot.set(i, baWhereNot[i]);
    if (level < parent->initialBaLevels())
       myBaWhereNot.coarsen(parent->refRatio(level));
    
    // we also have to restrict the creation of particles to non refined parts of the domain.
//    if (level == 0)
    Nyx::theDMPC()->InitCosmo1ppcMultiLevel(mf, dis_fac, vel_fac, 
		                  particleMass, 
				  part_dx, part_vx,
				  myBaWhereNot,
	                          level, parent->initialBaLevels()+1);
//    Nyx::theDMPC()->InitCosmo(mf, vel_fac, n_part, particleMass);
    //Nyx::theDMPC()->InitCosmo(mf, vel_fac, n_part, particleMass, part_dx, part_vx);

    MultiFab& S_new = get_level(level).get_new_data(State_Type);
    MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

#ifdef NUFLUID
    //copy density 
    S_new.copy(nuMf, nu_den, NuDens, 1);
    S_new.plus(1,            NuDens, 1, S_new.nGrow());
    S_new.mult(rhoNu,        NuDens, 1, S_new.nGrow());

    //copy velocities...
    S_new.copy(nuMf, nu_vx,  NuXMom, 3);

    //...and "transform" to momentum
    S_new.mult(vel_fac[0],   NuXMom, 1, S_new.nGrow());
    S_new.mult(vel_fac[1],   NuYMom, 1, S_new.nGrow());
    S_new.mult(vel_fac[2],   NuZMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuXMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuYMom, 1, S_new.nGrow());
    MultiFab::Multiply(S_new, S_new, NuDens, NuZMom, 1, S_new.nGrow());
#endif


    //initialize hydro
    //units for velocity should be okay
    if (do_hydro)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "Do hydro initialization..." << '\n';
        }

    	MultiFab& S_new = get_level(level).get_new_data(State_Type);
    	MultiFab& D_new = get_level(level).get_new_data(DiagEOS_Type);

        //Fill everything with old data...should only affect ghostzones, but
	//seems to have no effect...
	if (level > 0)
	{
           FillCoarsePatch(S_new, 0, 0,   State_Type, 0, S_new.nComp());
	   FillCoarsePatch(D_new, 0, 0, DiagEOS_Type, 0, D_new.nComp());
	}

     	//copy density 
     	S_new.setVal(0.);
     	S_new.copy(mf, baryon_den, Density, 1);
     	S_new.plus(1,     Density, 1, S_new.nGrow());
     	S_new.mult(rhoB,  Density, 1, S_new.nGrow());

//      //This block assigns "the same" density for the baryons as for the dm.
//      Vector<std::unique_ptr<MultiFab> > particle_mf;
//      Nyx::theDMPC()->AssignDensity(particle_mf);
//      particle_mf[0]->mult(comoving_OmB / comoving_OmD);
//      S_new.copy(*particle_mf[0], 0, Density, 1);

     	//copy velocities...
     	S_new.copy(mf, baryon_vx, Xmom, 3);

        //...and "transform" to momentum
     	S_new.mult(vel_fac[0], Xmom, 1, S_new.nGrow());
     	MultiFab::Multiply(S_new, S_new, Density, Xmom, 1, S_new.nGrow());
     	S_new.mult(vel_fac[1], Ymom, 1, S_new.nGrow());
     	MultiFab::Multiply(S_new, S_new, Density, Ymom, 1, S_new.nGrow());
     	S_new.mult(vel_fac[2], Zmom, 1, S_new.nGrow());
        MultiFab::Multiply(S_new, S_new, Density, Zmom, 1, S_new.nGrow());

        // Convert (rho X)_i to X_i before calling init_e_from_T
//      if (use_const_species == 0)
//          for (int i = 0; i < NumSpec; i++) 
//              MultiFab::Divide(S_new, S_new, Density, FirstSpec+i, 1, 0);

        Real tempInit = 0.021*(1.0+redshift)*(1.0+redshift);

        int ns = S_new.nComp();
        int nd = D_new.nComp();

        D_new.setVal(tempInit, Temp_comp);
        D_new.setVal(0.0, Ne_comp);
        if (inhomo_reion > 0)
            D_new.setVal(0.0, Zhi_comp);

#ifdef _OPENMP
#pragma omp parallel
#endif
     	for (MFIter mfi(S_new,true); mfi.isValid(); ++mfi)
     	{
     	    const Box& box = mfi.tilebox();
     	    const int* lo = box.loVect();
     	    const int* hi = box.hiVect();

     	    fort_init_e_from_t
      	       (BL_TO_FORTRAN(S_new[mfi]), &ns, 
      	        BL_TO_FORTRAN(D_new[mfi]), &nd, lo, hi, &old_a);
     	}     	

        // Convert X_i to (rho X)_i
        if (use_const_species == 0)
        {
           S_new.setVal(0.75, FirstSpec);
           S_new.setVal(0.25, FirstSpec+1);

           for (int i = 0; i < NumSpec; i++)
           {
              MultiFab::Multiply(S_new, S_new, Density, FirstSpec+i, 1, 0);
           }
        }

    }
    else
    {
       S_new.setVal(0.0, Density);
       D_new.setVal(-23, Temp_comp);
       D_new.setVal(-42, Ne_comp);
    }
       
    mf.clear();

#endif //ifdef GRAVITY
}
