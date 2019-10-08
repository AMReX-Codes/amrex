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
#include <MyParams.H>

void 
MyParams::test_parameters ()
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
