
#include <CNS.H>
#include <CNS_F.H>

using namespace amrex;

static Box the_same_box (const Box& b) { return b; }
//static Box grow_box_by_one (const Box& b) { return amrex::grow(b,1); }

using BndryFunc = StateDescriptor::BndryFunc;

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
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_ODD,  REFLECT_ODD,  REFLECT_ODD
};

static int tang_vel_bc[] =
{
    INT_DIR, EXT_DIR, FOEXTRAP, REFLECT_EVEN, REFLECT_EVEN, REFLECT_ODD
};

static
void
set_scalar_bc (BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    for (int i = 0; i < AMREX_SPACEDIM; i++)
    {
	bc.setLo(i,scalar_bc[lo_bc[i]]);
	bc.setHi(i,scalar_bc[hi_bc[i]]);
    }
}

static
void
set_x_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,norm_vel_bc[lo_bc[0]]);
    bc.setHi(0,norm_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_y_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)    
    bc.setLo(1,norm_vel_bc[lo_bc[1]]);
    bc.setHi(1,norm_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,tang_vel_bc[lo_bc[2]]);
    bc.setHi(2,tang_vel_bc[hi_bc[2]]);
#endif
}

static
void
set_z_vel_bc(BCRec& bc, const BCRec& phys_bc)
{
    const int* lo_bc = phys_bc.lo();
    const int* hi_bc = phys_bc.hi();
    bc.setLo(0,tang_vel_bc[lo_bc[0]]);
    bc.setHi(0,tang_vel_bc[hi_bc[0]]);
#if (AMREX_SPACEDIM >= 2)
    bc.setLo(1,tang_vel_bc[lo_bc[1]]);
    bc.setHi(1,tang_vel_bc[hi_bc[1]]);
#endif
#if (AMREX_SPACEDIM == 3)
    bc.setLo(2,norm_vel_bc[lo_bc[2]]);
    bc.setHi(2,norm_vel_bc[hi_bc[2]]);
#endif
}

void
CNS::variableSetUp ()
{
    read_params();

    cns_init_fort(phys_bc.lo(), phys_bc.hi(),
                  Interior,Inflow,Outflow,Symmetry,SlipWall,NoSlipWall,
                  Geometry::ProbLo(), Geometry::ProbHi());

    bool state_data_extrap = false;
    bool store_in_checkpoint = true;
    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
			   StateDescriptor::Point,NUM_GROW,NUM_STATE,
			   &eb_cell_cons_interp,state_data_extrap,store_in_checkpoint);

    Array<BCRec>       bcs(NUM_STATE);
    Array<std::string> name(NUM_STATE);
    BCRec bc;
    int cnt = 0;
    set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "density";
    cnt++; set_x_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "xmom";
    cnt++; set_y_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "ymom";
    cnt++; set_z_vel_bc(bc,phys_bc);  bcs[cnt] = bc; name[cnt] = "zmom";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_E";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "rho_e";
    cnt++; set_scalar_bc(bc,phys_bc); bcs[cnt] = bc; name[cnt] = "Temp";

    desc_lst.setComponent(State_Type,
			  Density,
			  name,
			  bcs,
			  BndryFunc(cns_denfill,cns_hypfill));

    StateDescriptor::setBndryFuncThreadSafety(true);

    // DEFINE DERIVED QUANTITIES

    // Pressure
    derive_lst.add("pressure",IndexType::TheCellType(),1,
                   cns_derpres,the_same_box);
    derive_lst.addComponent("pressure",desc_lst,State_Type,Eint,1);

    // Velocities
    derive_lst.add("x_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("x_velocity",desc_lst,State_Type,Xmom,1);

    derive_lst.add("y_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("y_velocity",desc_lst,State_Type,Ymom,1);

    derive_lst.add("z_velocity",IndexType::TheCellType(),1,
                   cns_dervel,the_same_box);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Density,1);
    derive_lst.addComponent("z_velocity",desc_lst,State_Type,Zmom,1);

    // Tagging 
    ErrorSetUp();
}

void
CNS::variableCleanUp ()
{
    desc_lst.clear();
}
