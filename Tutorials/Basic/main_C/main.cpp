//
// This example demonstrates:
//
// (1) By default, AMReX initializes MPI and uses MPI_COMM_WORLD as its communicator.
//     However, applications could choose to initialize MPI themselves and pass in an
//     existing communicator.
//
// (2) By default, AMReX treats command line arguments as inputs parameters.  The expected
//     format of argv is
//         executable inputs_file parm=value
//     Here, `excutable` is the filename of the executable, `inputs_file` is the file containing
//     runtime parameters used to build AMReX ParmParse database, and `parm=value` is an input
//     parameter that will override its value in `inputs_file`.  Both `inputs_file` and
//     `parm=value` are optional.  At most one `inputs_file` is allowed. Howeer, there can be
//     multiple `parm=value`s.
//
//     The parsing of the command line arguments is performed in amrex::Initialize.  Applications
//     can choose to skip command line parsing.  Applications can also provide a function that
//     adds parameters to AMReX ParmParse database.
//


#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

void add_parameters ();
void test_parameters ();

int main(int argc, char* argv[])
{
    // We choose to call MPI_Init ourselves, instead of letting AMRex do it.
    // Then we are going to pass a communicatior, MPI_COMM_WORLD, (or another) to AMReX.
    MPI_Init(&argc, &argv);
    
    bool build_parm_parse_from_command_line = false;
    amrex::Initialize(argc, argv, build_parm_parse_from_command_line,
                      MPI_COMM_WORLD, add_parameters);  // last three arguments are optional

    test_parameters();

    // ...

    amrex::Finalize();

    MPI_Finalize();  // We have to call this because we called MPI_Init.
}

void add_parameters ()
{
    {  // prefix "amrex"
        amrex::ParmParse pp("amrex");
        pp.add("fpe_trap_invalid",1); //  turn on NaN trapping, which is off by default.
    }

    {  // anonymous prefix
        amrex::ParmParse pp;
        pp.add("an_int_scalar", 2);            // integer scalar: an_int_scalar
        pp.add("a_bool_scalar", true);         // logical scalar: a_bool_scalar
        pp.addarr("a_real_array",              // real array: a_real_array
                  std::vector<amrex::Real>{1.,2.,3.}); 
    }

    {
        // prefix "a_prefix"
        amrex::ParmParse pp("a_prefix");
        pp.addarr("an_int_array",              // integer array: a_prefix.an_int_array
                  std::vector<int>{2, 3, 4});      
        amrex::Real x = 3.14;
        pp.add("a_real_scalar", x);            // real scalar  : a_prefix.a_real_scalar
        pp.add("a_string", std::string{"vonNeumann"});  // string: a_prefix.a_string
    }
}

void test_parameters ()
{
    {
        amrex::ParmParse pp;
        int i;
        bool b;
        std::vector<amrex::Real> ra;
        pp.get("an_int_scalar", i);
        pp.get("a_bool_scalar",b);
        pp.getarr("a_real_array", ra);
        amrex::Print() << "an_int_scalar = " << i << "\n"
                       << "a_bool_scalar = " << b << "\n";
        amrex::Print() << "a_real_array = ";
        for (auto x : ra) {
            amrex::Print() << x << " ";
        }
        amrex::Print() << "\n";
    }

    {
        amrex::ParmParse pp("a_prefix");
        std::vector<int> ia;
        amrex::Real r;
        std::string s;
        pp.getarr("an_int_array", ia);
        pp.get("a_real_scalar", r);
        pp.get("a_string", s);
        amrex::Print() << "an_int_array = ";
        for (auto x : ia) {
            amrex::Print() << x << " ";
        }
        amrex::Print() << "\n";
        amrex::Print() << "a_prefix.a_real_scalar = " << r << "\n"
                       << "a_prefix.a_string = " << s << "\n";
    }
}
