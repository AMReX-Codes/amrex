
#include <Adv.H>
#include <Adv_F.H>
#include <ParmParse.H>

void
Adv::variableSetUp ()
{
    BL_ASSERT(desc_lst.size() == 0);

    // Get options, set phys_bc
    read_params();

    desc_lst.addDescriptor(State_Type,IndexType::TheCellType(),
                           StateDescriptor::Point,0,NUM_STATE,
			   &cell_cons_interp);

    desc_lst.setComponent(State_Type, 0, "phi", BCRec(), StateDescriptor::BndryFunc());

    //
    // read taggin parameters from probin file
    //

    std::string probin_file("probin");

    ParmParse ppa("amr");
    ppa.query("probin_file",probin_file);

    int probin_file_length = probin_file.length();
    Array<int> probin_file_name(probin_file_length);

    for (int i = 0; i < probin_file_length; i++)
	probin_file_name[i] = probin_file[i];

     BL_FORT_PROC_CALL(GET_TAGGING_PARAMS, get_tagging_params)
	(probin_file_name.dataPtr(), &probin_file_length);
}
