#include <AMReX.H>
#include <AMReX_ParmParse.H>

#include <exception>
#include <string>
#include <vector>

namespace {

constexpr int         ref_int    = 42;
constexpr amrex::Real ref_real   = 2.5;
constexpr int         ref_qadd   = 7;
constexpr int         ref_file   = 1234;
constexpr int         ref_pfac   = 11;

std::string ctx (char const* what, int cycle)
{
    return std::string(what) + " leaked into cycle " + std::to_string(cycle)
         + ": amrex::Finalize() did not clear the ParmParse table";
}

/** Everything an earlier cycle put into the ParmParse table must be gone.
 *
 * amrex::Finalize() clears the table, so a new Initialize() starts from a
 * clean slate no matter how AMReX was initialized. Without this, a process
 * that runs several simulations back to back silently inherits the previous
 * run's parameters.
 */
void assert_table_is_clean (int cycle)
{
    {
        amrex::ParmParse pp("reinit");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("an_int"),    ctx("reinit.an_int", cycle));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("a_real"),    ctx("reinit.a_real", cycle));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("a_string"),  ctx("reinit.a_string", cycle));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("an_array"),  ctx("reinit.an_array", cycle));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("queried"),   ctx("reinit.queried", cycle));
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("expr"),      ctx("reinit.expr", cycle));
    }
    {
        // entries pulled in with ParmParse::addfile must not survive either
        amrex::ParmParse pp("reinit_file");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("from_file"), ctx("reinit_file.from_file", cycle));
    }
    {
        amrex::ParmParse pp("reinit_parser");
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!pp.contains("pfac"),      ctx("reinit_parser.pfac", cycle));
    }
}

/** Entries from the inputs file must come back in every cycle.
 *
 * Counterpart to assert_table_is_clean(): the table is cleared on Finalize,
 * but Initialize re-reads the command line and the inputs file, so anything
 * the user passed there is present again.
 */
void assert_inputs_file_entries_present (int cycle)
{
    amrex::ParmParse pp("amrex");
    int verbose = -1;
    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(pp.query("verbose", verbose) && verbose == 2,
        "amrex.verbose from the inputs file is missing in cycle " + std::to_string(cycle));
}

/** Add table entries the way a running application does. */
void add_runtime_entries ()
{
    {
        amrex::ParmParse pp("reinit");
        pp.add("an_int", ref_int);
        pp.add("a_real", ref_real);
        pp.add("a_string", std::string("hello"));
        pp.addarr("an_array", std::vector<int>{1,2,3});

        // queryAdd writes its default into the table when the key is absent
        int queried = ref_qadd;
        pp.queryAdd("queried", queried);
        AMREX_ALWAYS_ASSERT(queried == ref_qadd);
    }

    // read them back within the same cycle
    {
        amrex::ParmParse pp("reinit");
        int an_int = 0;
        amrex::Real a_real = 0.0;
        std::string a_string;
        std::vector<int> an_array;
        AMREX_ALWAYS_ASSERT(pp.query("an_int", an_int)     && an_int == ref_int);
        AMREX_ALWAYS_ASSERT(pp.query("a_real", a_real)     && a_real == ref_real);
        AMREX_ALWAYS_ASSERT(pp.query("a_string", a_string) && a_string == "hello");
        AMREX_ALWAYS_ASSERT(pp.queryarr("an_array", an_array) && an_array.size() == 3);
    }

    amrex::ParmParse::addfile("inputs.addfile");
    {
        amrex::ParmParse pp("reinit_file");
        int from_file = 0;
        AMREX_ALWAYS_ASSERT(pp.query("from_file", from_file) && from_file == ref_file);
    }
}

/** The parser prefix set with SetParserPrefix must not survive a cycle.
 *
 * It is file-scope state next to the table, so Finalize has to reset it too.
 * Applications set it from the func_parm_parse hook, i.e. once per cycle.
 */
void check_parser_prefix (bool expect_prefix_set, int cycle)
{
    // only exercised when errors throw; otherwise the negative case would abort
    if (!amrex::system::throw_exception) { return; }

    amrex::ParmParse pp;
    pp.add("reinit_parser.pfac", ref_pfac);
    pp.add("reinit.expr", std::string("2*pfac"));

    int value = 0;
    bool resolved = false;
    try {
        resolved = pp.queryWithParser("reinit.expr", value);
    } catch (std::exception const&) {
        resolved = false;  // unknown symbol: no parser prefix in effect
    }

    if (expect_prefix_set) {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(resolved && value == 2*ref_pfac,
            "parser prefix set in this cycle did not resolve the symbol");
    } else {
        AMREX_ALWAYS_ASSERT_WITH_MESSAGE(!resolved,
            "ParmParse::ParserPrefix leaked into cycle " + std::to_string(cycle));
    }
}

}

int main (int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    int cycle = 0;

    // Cycles driven by the command line (inputs file): argc > 1.
    for (int i = 0; i < 3; ++i)
    {
        ++cycle;
        amrex::Initialize(argc,argv);

        assert_table_is_clean(cycle);
        assert_inputs_file_entries_present(cycle);

        // the first cycle sets the parser prefix, later ones must not see it
        bool const set_prefix = (i == 0);
        if (set_prefix) { amrex::ParmParse::SetParserPrefix("reinit_parser"); }
        check_parser_prefix(set_prefix, cycle);

        add_runtime_entries();
        amrex::Finalize();
    }

    // Cycles without a command line: argc == 0 with build_parm_parse = true is
    // how AMReX is initialized from a library or from language bindings.
    for (int i = 0; i < 2; ++i)
    {
        ++cycle;
        int argc0 = 0;
        char** argv0 = nullptr;
        amrex::Initialize(argc0,argv0,true);

        assert_table_is_clean(cycle);

        add_runtime_entries();
        amrex::Finalize();
    }

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
