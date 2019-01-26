
#include "AMReX_LevelBld.H"

#include "Nyx.H"
#include "Nyx_F.H"
#include "Derive_F.H"

using namespace amrex;
using std::string;

static Box the_same_box(const Box& b)
{
    return b;
}

static Box grow_box_by_one(const Box& b)
{
    return amrex::grow(b, 1);
}

typedef StateDescriptor::BndryFunc BndryFunc;

//
// Components are:
//  Interior, Inflow, Outflow,  Symmetry,     SlipWall,     NoSlipWall
//
static int scalar_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static int norm_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD, REFLECT_ODD, REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_EVEN
};

static
void
set_scalar_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        bc.setLo(i, scalar_bc[lo_bc[i]]);
        bc.setHi(i, scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, norm_vel_bc[lo_bc[0]]);
    bc.setHi(0, norm_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, norm_vel_bc[lo_bc[1]]);
    bc.setHi(1, norm_vel_bc[hi_bc[1]]);
    bc.setLo(2, tang_vel_bc[lo_bc[2]]);
    bc.setHi(2, tang_vel_bc[hi_bc[2]]);
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0, tang_vel_bc[lo_bc[0]]);
    bc.setHi(0, tang_vel_bc[hi_bc[0]]);
    bc.setLo(1, tang_vel_bc[lo_bc[1]]);
    bc.setHi(1, tang_vel_bc[hi_bc[1]]);
    bc.setLo(2, norm_vel_bc[lo_bc[2]]);
    bc.setHi(2, norm_vel_bc[hi_bc[2]]);
}

void
Nyx::variable_setup()
{

  // initialize the start time for our CPU-time tracker
  startCPUTime = ParallelDescriptor::second();

    BL_ASSERT(desc_lst.size() == 0);

    // Initialize the network
    network_init();

    // Get options, set phys_bc
    read_params();

#ifdef NO_HYDRO
    no_hydro_setup();

#else
    if (do_hydro == 1)
    {
       hydro_setup();
    }
#ifdef GRAVITY
    else
    {
       no_hydro_setup();
    }
#endif
#endif

    //
    // DEFINE ERROR ESTIMATION QUANTITIES
    //
    error_setup();
}

#ifndef NO_HYDRO
void
Nyx::hydro_setup()
{
    //
    // Set number of state variables and pointers to components
    //
    Density = 0;
    Xmom    = 1;
    Ymom    = 2;
    Zmom    = 3;
    Eden    = 4;
    Eint    = 5;
    int cnt = 6;

    NumAdv = 0;
    if (NumAdv > 0)
    {
        FirstAdv = cnt;
        cnt += NumAdv;
    }

    int NDIAG_C;
    Temp_comp = 0;
      Ne_comp = 1;
    if (inhomo_reion > 0)
    {
        NDIAG_C  = 3;
        Zhi_comp = 2;
    } else {
        NDIAG_C  = 2;
    }

    int dm = BL_SPACEDIM;

    // Get the number of species from the network model.
    fort_get_num_spec(&NumSpec);

    if (use_const_species == 0)
    {
        if (NumSpec > 0)
        {
            FirstSpec = cnt;
            cnt += NumSpec;
        }

        // Get the number of auxiliary quantities from the network model.
        fort_get_num_aux(&NumAux);

        if (NumAux > 0)
        {
            FirstAux = cnt;
            cnt += NumAux;
        }
    }

    NUM_STATE = cnt;

    // Define NUM_GROW from the f90 module.
    fort_get_method_params(&NUM_GROW);

    // Note that we must set NDIAG_C before we call set_method_params because
    // we use the C++ value to set the Fortran value
    fort_set_method_params
        (dm, NumAdv, NDIAG_C, do_hydro, ppm_type, ppm_reference,
         ppm_flatten_before_integrals,
         use_colglaz, use_flattening, corner_coupling, version_2,
         use_const_species, gamma, normalize_species,
         heat_cool_type, inhomo_reion);

#ifdef HEATCOOL
    fort_tabulate_rates();
#endif

    if (use_const_species == 1)
        fort_set_eos_params(h_species, he_species);

    int coord_type = Geometry::Coord();
    fort_set_problem_params
         (dm, phys_bc.lo(), phys_bc.hi(), Outflow, Symmetry, coord_type);

    Interpolater* interp = &cell_cons_interp;

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, NUM_STATE, interp,
                           state_data_extrap, store_in_checkpoint);

    // This has two components: Temperature and Ne
    desc_lst.addDescriptor(DiagEOS_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, NDIAG_C, interp,
                           state_data_extrap, store_in_checkpoint);

#ifdef SDC
    // This only has one component -- the update to rho_e from reactions
    store_in_checkpoint = true;
    desc_lst.addDescriptor(SDC_IR_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, interp,
                           state_data_extrap, store_in_checkpoint);
#endif

#ifdef GRAVITY
    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, BL_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);
#endif

    Vector<BCRec> bcs(NUM_STATE);
    Vector<std::string> name(NUM_STATE);

    BCRec bc;
    cnt = 0;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "density";
    cnt++;
    set_x_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "xmom";
    cnt++;
    set_y_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "ymom";
    cnt++;
    set_z_vel_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "zmom";
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "rho_E";
    cnt++;
    set_scalar_bc(bc, phys_bc);  bcs[cnt] = bc;  name[cnt] = "rho_e";

    for (int i = 0; i < NumAdv; ++i)
    {
        cnt++;
        set_scalar_bc(bc, phys_bc);
        bcs[cnt]  = bc;
        name[cnt] = amrex::Concatenate("adv_", i, 1);
    }

    // Get the species names from the network model.
    Vector<std::string> spec_names(NumSpec);

    for (int i = 0; i < NumSpec; i++)
    {
        int len = 20;
        Vector<int> int_spec_names(len);

        // This call return the actual length of each string in "len"
        fort_get_spec_names
            (int_spec_names.dataPtr(), &i, &len);

        for (int j = 0; j < len; j++)
            spec_names[i].push_back(int_spec_names[j]);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumSpec << " Species: ";
        for (int i = 0; i < NumSpec; i++)
           std::cout << spec_names[i] << ' ' << ' ';
        std::cout << '\n';
    }

    if (use_const_species == 0)
    {
        for (int i = 0; i < NumSpec; ++i)
        {
            cnt++;
            set_scalar_bc(bc,phys_bc);
            bcs[cnt] = bc;
                name[cnt] = "rho_" + spec_names[i];
         }
    }

    // Get the auxiliary names from the network model.
    Vector<std::string> aux_names(NumAux);

    for (int i = 0; i < NumAux; i++)
    {
        int len = 20;
        Vector<int> int_aux_names(len);

        // This call return the actual length of each string in "len"
        fort_get_aux_names
            (int_aux_names.dataPtr(), &i, &len);

        for (int j = 0; j < len; j++)
            aux_names[i].push_back(int_aux_names[j]);
    }

    if (ParallelDescriptor::IOProcessor())
    {
        std::cout << NumAux << " Auxiliary Variables: ";
        for (int i = 0; i < NumAux; i++)
           std::cout << aux_names[i] << ' ' << ' ';
        std::cout << '\n';
    }

    if (use_const_species == 0)
    {
        for (int i = 0; i < NumAux; ++i)
        {
            cnt++;
            set_scalar_bc(bc, phys_bc);
            bcs[cnt] = bc;
            name[cnt] = "rho_" + aux_names[i];
        }
    }

    desc_lst.setComponent(State_Type, Density, name, bcs,
                          BndryFunc(denfill,hypfill));

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(DiagEOS_Type, 0, "Temp", bc,
                          BndryFunc(generic_fill));
    desc_lst.setComponent(DiagEOS_Type, 1, "Ne", bc,
                          BndryFunc(generic_fill));

    if (inhomo_reion > 0) {
       desc_lst.setComponent(DiagEOS_Type, 2, "Z_HI", bc,
                             BndryFunc(generic_fill));
    }

#ifdef SDC
    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(SDC_IR_Type, 0, "I_R", bc,
                          BndryFunc(generic_fill));
#endif

#ifdef GRAVITY
    if (do_grav)
    {
        set_scalar_bc(bc, phys_bc);
        desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                              BndryFunc(generic_fill));
        set_x_vel_bc(bc, phys_bc);
        desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                              BndryFunc(generic_fill));
       set_y_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                             BndryFunc(generic_fill));
       set_z_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                             BndryFunc(generic_fill));
    }
#endif

    //
    // DEFINE DERIVED QUANTITIES
    //
    // Pressure
    //
    derive_lst.add("pressure", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERPRES, derpres), the_same_box);
    derive_lst.addComponent("pressure", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Kinetic energy
    //
    derive_lst.add("kineng", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERKINENG, derkineng), the_same_box);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("kineng", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Sound speed (c)
    //
    derive_lst.add("soundspeed", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERSOUNDSPEED, dersoundspeed),
                   the_same_box);
    derive_lst.addComponent("soundspeed", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Mach number(M)
    //
    derive_lst.add("MachNumber", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERMACHNUMBER, dermachnumber),
                   the_same_box);
    derive_lst.addComponent("MachNumber", desc_lst, State_Type, Density,
                            NUM_STATE);

    //
    // Gravitational forcing
    //
#ifdef GRAVITY
    //derive_lst.add("rhog",IndexType::TheCellType(),1,
    //               BL_FORT_PROC_CALL(RHOG,rhog),the_same_box);
    //derive_lst.addComponent("rhog",desc_lst,State_Type,Density,1);
    //derive_lst.addComponent("rhog",desc_lst,Gravity_Type,0,BL_SPACEDIM);
#endif

    //
    // Div(u)
    //
    derive_lst.add("divu", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERDIVU, derdivu), grow_box_by_one);
    derive_lst.addComponent("divu", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("divu", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Internal energy as derived from rho*E, part of the state
    //
    derive_lst.add("eint_E", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DEREINT1, dereint1), the_same_box);
    derive_lst.addComponent("eint_E", desc_lst, State_Type, Density, NUM_STATE);

    //
    // Internal energy as derived from rho*e, part of the state
    //
    derive_lst.add("eint_e", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DEREINT2, dereint2), the_same_box);
    derive_lst.addComponent("eint_e", desc_lst, State_Type, Density, NUM_STATE);

    //
    // Log(density)
    //
    derive_lst.add("logden", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERLOGDEN, derlogden), the_same_box);
    derive_lst.addComponent("logden", desc_lst, State_Type, Density, 1);

    derive_lst.add("StateErr", IndexType::TheCellType(), 3,
                   BL_FORT_PROC_CALL(DERSTATE, derstate),
                   grow_box_by_one);
    derive_lst.addComponent("StateErr", desc_lst,   State_Type, Density, 1);
    derive_lst.addComponent("StateErr", desc_lst, DiagEOS_Type, Temp_comp, 1);
    if (use_const_species == 0)
        derive_lst.addComponent("StateErr", desc_lst, State_Type, FirstSpec, 1);

    //
    // X from rhoX
    //
    if (use_const_species == 0)
    {
        for (int i = 0; i < NumSpec; i++)
        {
            string spec_string = "X(" + spec_names[i] + ")";

            derive_lst.add(spec_string, IndexType::TheCellType(), 1,
                           BL_FORT_PROC_CALL(DERSPEC, derspec), the_same_box);
            derive_lst.addComponent(spec_string, desc_lst, State_Type, Density, 1);
            derive_lst.addComponent(spec_string, desc_lst, State_Type,
                                    FirstSpec + i, 1);
        }
    }

    //
    // Forcing
    //
    derive_lst.add("forcex", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERFORCEX, derforcex), the_same_box);
    derive_lst.addComponent("forcex", desc_lst, State_Type, Density, 1);

    derive_lst.add("forcey", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERFORCEY, derforcey), the_same_box);
    derive_lst.addComponent("forcey", desc_lst, State_Type, Density, 1);

    derive_lst.add("forcez", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERFORCEZ, derforcez), the_same_box);
    derive_lst.addComponent("forcez", desc_lst, State_Type, Density, 1);

    //
    // Velocities
    //
    derive_lst.add("x_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERVEL, dervel), the_same_box);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("x_velocity", desc_lst, State_Type, Xmom, 1);

    derive_lst.add("y_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERVEL, dervel), the_same_box);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("y_velocity", desc_lst, State_Type, Ymom, 1);

    derive_lst.add("z_velocity", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERVEL, dervel), the_same_box);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("z_velocity", desc_lst, State_Type, Zmom, 1);

    //
    // Magnitude of velocity.
    //
    derive_lst.add("magvel", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERMAGVEL, dermagvel), the_same_box);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Density, 1);
    derive_lst.addComponent("magvel", desc_lst, State_Type, Xmom, BL_SPACEDIM);

    //
    // Magnitude of vorticity.
    //
    derive_lst.add("magvort",IndexType::TheCellType(),1,
                   BL_FORT_PROC_CALL(DERMAGVORT,dermagvort),grow_box_by_one);
    // Here we exploit the fact that Xmom = Density + 1
    //   in order to use the correct interpolation.
    if (Xmom != Density+1)
       amrex::Error("We are assuming Xmom = Density + 1 in Nyx_setup.cpp");
    derive_lst.addComponent("magvort",desc_lst,State_Type,Density,BL_SPACEDIM+1);

    //
    // Magnitude of momentum.
    //
    derive_lst.add("magmom", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERMAGMOM, dermagmom), the_same_box);
    derive_lst.addComponent("magmom", desc_lst, State_Type, Xmom, BL_SPACEDIM);

#ifdef GRAVITY
    derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERMAGGRAV, dermaggrav),
                   the_same_box);
    derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, BL_SPACEDIM);
#endif

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, State_Type, Density, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, State_Type,
                            Density, 1);

#ifdef AGN
    derive_lst.add("agn_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("agn_particle_count", desc_lst, State_Type, Density, 1);

    derive_lst.add("agn_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("agn_mass_density", desc_lst, State_Type,
                            Density, 1);
#endif

#ifdef NEUTRINO_PARTICLES
    derive_lst.add("neutrino_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("neutrino_particle_count", desc_lst, State_Type, Density, 1);

    derive_lst.add("neutrino_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("neutrino_mass_density", desc_lst, State_Type,
                            Density, 1);
#endif

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, State_Type,
                            Density, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, State_Type,
                            Density, 1);
    derive_lst.add("Rank", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);

    if (use_const_species == 0)
    {
        for (int i = 0; i < NumSpec; i++)
        {
            derive_lst.add(spec_names[i], IndexType::TheCellType(), 1,
                           BL_FORT_PROC_CALL(DERSPEC, derspec), the_same_box);
            derive_lst.addComponent(spec_names[i], desc_lst, State_Type, Density, 1);
            derive_lst.addComponent(spec_names[i], desc_lst, State_Type,
                                    FirstSpec + i, 1);
        }

        for (int i = 0; i < NumAux; i++)
        {
            derive_lst.add(aux_names[i], IndexType::TheCellType(), 1,
                           BL_FORT_PROC_CALL(DERSPEC, derspec), the_same_box);
            derive_lst.addComponent(aux_names[i], desc_lst, State_Type, Density, 1);
            derive_lst.addComponent(aux_names[i], desc_lst, State_Type, FirstAux+i, 1);
        }
    }
}
#endif

#ifdef GRAVITY
void
Nyx::no_hydro_setup()
{
    int dm = BL_SPACEDIM;

    Density = 0;
    NUM_STATE = 1;

    int NDIAG_C = -1;

    // Define NUM_GROW from the f90 module.
    fort_get_method_params(&NUM_GROW);

    fort_set_method_params
        (dm, NumAdv, NDIAG_C, do_hydro, ppm_type, ppm_reference,
         ppm_flatten_before_integrals,
         use_colglaz, use_flattening, corner_coupling, version_2,
         use_const_species, gamma, normalize_species,
         heat_cool_type, inhomo_reion);

    fort_tabulate_rates();

    int coord_type = Geometry::Coord();
    fort_set_problem_params(dm, phys_bc.lo(), phys_bc.hi(), Outflow, Symmetry, coord_type);

    // Note that the default is state_data_extrap = false,
    // store_in_checkpoint = true.  We only need to put these in
    // explicitly if we want to do something different,
    // like not store the state data in a checkpoint directory
    bool state_data_extrap = false;
    bool store_in_checkpoint;

    BCRec bc;

#ifndef NO_HYDRO
    // We have to create these anyway because the StateTypes are defined at compile time.
    // However, we define them with only one component each because we don't actually use them.
    store_in_checkpoint = false;
    // This has only one dummy components
    desc_lst.addDescriptor(State_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, &cell_cons_interp,
                           state_data_extrap, store_in_checkpoint);

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(State_Type, 0, "density", bc,
                          BndryFunc(generic_fill));

    // This has only one dummy components
    store_in_checkpoint = false;
    desc_lst.addDescriptor(DiagEOS_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1, &cell_cons_interp,
                           state_data_extrap, store_in_checkpoint);

    set_scalar_bc(bc, phys_bc);
    desc_lst.setComponent(DiagEOS_Type, 0, "Temp", bc,
                          BndryFunc(generic_fill));
#endif

    store_in_checkpoint = true;
    desc_lst.addDescriptor(PhiGrav_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, 1,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    store_in_checkpoint = false;
    desc_lst.addDescriptor(Gravity_Type, IndexType::TheCellType(),
                           StateDescriptor::Point, 1, BL_SPACEDIM,
                           &cell_cons_interp, state_data_extrap,
                           store_in_checkpoint);

    if (do_grav)
    {
       set_scalar_bc(bc, phys_bc);
       desc_lst.setComponent(PhiGrav_Type, 0, "phi_grav", bc,
                             BndryFunc(generic_fill));
       set_x_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 0, "grav_x", bc,
                             BndryFunc(generic_fill));
       set_y_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 1, "grav_y", bc,
                             BndryFunc(generic_fill));
       set_z_vel_bc(bc, phys_bc);
       desc_lst.setComponent(Gravity_Type, 2, "grav_z", bc,
                             BndryFunc(generic_fill));

       derive_lst.add("maggrav", IndexType::TheCellType(), 1,
                      BL_FORT_PROC_CALL(DERMAGGRAV, dermaggrav),
                      the_same_box);
       derive_lst.addComponent("maggrav", desc_lst, Gravity_Type, 0, BL_SPACEDIM);
    }

    //
    // We want a derived type that corresponds to the number of particles in
    // each cell. We only intend to use it in plotfiles for debugging purposes.
    // We'll just use the DERNULL since don't do anything in fortran for now.
    // We'll actually set the values in `writePlotFile()`.
    //
    derive_lst.add("particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("particle_count", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("particle_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("particle_mass_density", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("total_particle_count", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), the_same_box);
    derive_lst.addComponent("total_particle_count", desc_lst, PhiGrav_Type, 0, 1);

    derive_lst.add("total_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("total_density", desc_lst, PhiGrav_Type, 0, 1);

#ifdef AGN
    derive_lst.add("agn_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("agn_mass_density", desc_lst, Gravity_Type, 0, 1);
#endif
#ifdef NEUTRINO_PARTICLES
    derive_lst.add("neutrino_mass_density", IndexType::TheCellType(), 1,
                   BL_FORT_PROC_CALL(DERNULL, dernull), grow_box_by_one);
    derive_lst.addComponent("neutrino_mass_density", desc_lst, Gravity_Type, 0, 1);
#endif
}
#endif

#ifdef AMREX_USE_CVODE
void
Nyx::set_simd_width(const int simd_width)
{
    set_simd(&simd_width);
}

void
Nyx::alloc_simd_vec()
{
    fort_alloc_simd_vec();
}

void
Nyx::dealloc_simd_vec()
{
    fort_dealloc_simd_vec();
}
#endif
