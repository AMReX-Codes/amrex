#include <AMReX.H>
#include <AMReX_Box.H>
#include <AMReX_Utility.H>
#include <AMReX_Print.H>
#include <AMReX_ParmParse.H>

using namespace amrex;

int main(int argc, char* argv[])
{
#if !defined(_WIN32)
    if (! std::getenv("AMREX_DEFAULT_INIT")) {
        setenv("AMREX_DEFAULT_INIT",
               R"(amrex.envfoo=0 amrex.envbar=1 amrex.envabc=1 2 3 amrex.envstr="a b c")", 1);
    }
#endif

    amrex::Initialize(argc,argv);
    {
        ParmParse::SetParserPrefix("physical_constants");
        ParmParse pp("physical_constants");
        pp.add("c", 299792458.);
        pp.add("pi", 3.14159265358979323846);
    }
    {
        ParmParse pp;
        int val;
        pp.query("dAx_x/dx(x,y,t,zeval)", val);
        AMREX_ALWAYS_ASSERT(val == 12);
    }
    {
        ParmParse pp;

        std::string name;
        pp.query("name", name);
        AMREX_ALWAYS_ASSERT(name == "I am w");
        pp.query("name", name, 1);
        AMREX_ALWAYS_ASSERT(name == "line 2");

        Box box;
        pp.query("b", box);
        AMREX_ALWAYS_ASSERT(box == Box(IntVect(AMREX_D_DECL(1,2,3)),
                                       IntVect(AMREX_D_DECL(7,8,9)),
                                       IntVect(AMREX_D_DECL(1,0,1))));

        double f0 = -1;
        pp.query("f", f0);
        AMREX_ALWAYS_ASSERT(f0 == 7);

        std::vector<int> f;
        pp.queryarr("f", f);
        AMREX_ALWAYS_ASSERT(f[0] == 7 && f[1] == 99 && f[2] == 11);

        std::vector<double> g;
        pp.queryarr("g", g);
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(g[0], 7.2) &&
                            amrex::almostEqual(g[1], 11.6));

        double w;
        pp.query("w", w);
        AMREX_ALWAYS_ASSERT(w == 1);
        pp.queryWithParser("w", w);
        AMREX_ALWAYS_ASSERT(w == -1);
    }
    {
        ParmParse pp("amrex", "my_constants");
        double foo = -1, bar = -2, bar2 = -3;
        pp.getWithParser("foo", foo);
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(foo, 6.0-std::sqrt(299792458.)));
        pp.get("bar", bar);
        AMREX_ALWAYS_ASSERT(foo == bar);
        pp.get("bar2", bar2);
        AMREX_ALWAYS_ASSERT(bar == bar2);
    }
    {
        ParmParse pp;
        std::array<double,3> prob_lo, prob_hi;
        pp.get("geom.prob_lo", prob_lo);
        pp.get("geom.prob_hi", prob_hi);
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(prob_lo[0], -1.0) &&
                            amrex::almostEqual(prob_lo[1], -1.0) &&
                            amrex::almostEqual(prob_lo[2], -1.0) &&
                            amrex::almostEqual(prob_hi[0],  1.0) &&
                            amrex::almostEqual(prob_hi[1],  1.0) &&
                            amrex::almostEqual(prob_hi[2],  1.0));
    }
    {
        ParmParse pp;
        auto parser = pp.makeParser("pi*x+c*y", {"x","y"});
        auto exe = parser.compile<2>();
        AMREX_ALWAYS_ASSERT(amrex::almostEqual(3.14159265358979323846+299792458.,
                                               exe(1.0,1.0)) &&
                            amrex::almostEqual(3.14159265358979323846, exe(1.0,0.0)) &&
                            amrex::almostEqual(299792458., exe(0.0, 1.0)));
    }
    {
        ParmParse pp;
        long long int i = 123456789012345;
        long long int j = 0;
        pp.get("long_int_1", j);
        AMREX_ALWAYS_ASSERT(i==j);
        pp.get("long_int_2", j);
        AMREX_ALWAYS_ASSERT(i==j);
        pp.get("long_int_3", j);
        AMREX_ALWAYS_ASSERT(i==j);
    }
    try
    {
        ParmParse pp("code");
        int a = 0;
        pp.query("a",a);
        amrex::Abort("Should not get here, because query should raise an exception");
    } catch (std::runtime_error const& e) {
        // Runtime error as expected
        amrex::ignore_unused(e);
    }
    {
        int max_steps = -1;
        ParmParse pp("", "my_constants");
        pp.query("max_steps", max_steps);
        AMREX_ALWAYS_ASSERT(max_steps == 40);
        int warpx_max_steps = -1;
        pp.query("warpx.max_steps", warpx_max_steps);
        AMREX_ALWAYS_ASSERT(max_steps == 40);
    }
    {
        ParmParse::SetParserPrefix("my_constants");
        ParmParse pp;

        int ny = 0;
        pp.queryAsDouble("ny", ny);
        AMREX_ALWAYS_ASSERT(ny == 64);

        Array<int,3> n_cell{0,0,0};
        pp.queryarrAsDouble("n_cell", 3, n_cell.data());
        AMREX_ALWAYS_ASSERT(n_cell[0] == 64 && n_cell[1] == 64 && n_cell[2] == 64);
    }
    {
        ParmParse pp;
        bool b_do_this = false;
        pp.queryAsDouble("do_this", b_do_this);
        AMREX_ALWAYS_ASSERT(b_do_this);

        std::optional<int> o_do_this;
        pp.queryAsDouble("do_this", o_do_this);
        AMREX_ALWAYS_ASSERT(o_do_this.has_value() && o_do_this.value());

        std::optional<int> o_do_that;
        pp.queryAsDouble("do_that", o_do_that);
        AMREX_ALWAYS_ASSERT(!o_do_that.has_value());
    }
    {
        ParmParse pp;
        bool my_bool_flag_1 = false;
        bool my_bool_flag_2 = false;
        pp.queryAddWithParser("my_bool_flag", my_bool_flag_1);
        pp.query("my_bool_flag", my_bool_flag_2);
        AMREX_ALWAYS_ASSERT(my_bool_flag_1 && my_bool_flag_2);
    }
    {
        ParmParse pp;
        std::string line;
        pp.queryline("my_string_line", line);
        AMREX_ALWAYS_ASSERT(line == "a b c");
        line.clear();
        pp.getline("my_string_line", line);
        AMREX_ALWAYS_ASSERT(line == "a b c");
    }
#if !defined(_WIN32)
    {
        int envfoo, envbar;
        std::vector<int> envabc;
        std::string envstr;
        ParmParse pp("amrex");
        pp.get("envfoo", envfoo);
        pp.get("envbar", envbar);
        pp.getarr("envabc", envabc);
        pp.get("envstr", envstr);
        AMREX_ALWAYS_ASSERT(envfoo == 0 && envbar == 1 &&
                            envabc.size() == 3 &&
                            envabc[0] == 1 && envabc[1] == 2 && envabc[2] == 3 &&
                            envstr == "a b c");
    }
#endif
    {
        ParmParse pp("t");
        std::vector<std::vector<double>> table;
        pp.querytable("table", table);
        std::vector<std::vector<int>> table2;
        pp.gettable("table2", table2);
        AMREX_ALWAYS_ASSERT(table.size() == 4 && table2.size() == 4);
        for (int irow = 0; irow < 4; ++irow) {
            AMREX_ALWAYS_ASSERT(table[irow].size() == 3 && table2[irow].size() == 3);
            for (int icol = 0; icol < 3; ++icol) {
                AMREX_ALWAYS_ASSERT(table [irow][icol] == (irow+1)*10.+icol+1 &&
                                    table2[irow][icol] == (irow+1)*10 +icol+1);
            }
        }
    }
    {
        amrex::Print() << "SUCCESS\n";
    }
    amrex::Finalize();
}
