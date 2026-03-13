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

        std::vector<std::string> sa;
        std::vector<std::string> sb;
        pp.getarr("sa", sa);
        pp.getarr("sa", sb);
        AMREX_ALWAYS_ASSERT(sa == sb && (sa == std::vector<std::string>{"abc","xyz","123"}));

        Box box;
        pp.query("b", box);
        AMREX_ALWAYS_ASSERT(box == Box(IntVect(AMREX_D_DECL(1,2,3)),
                                       IntVect(AMREX_D_DECL(7,8,9)),
                                       IntVect(AMREX_D_DECL(1,0,1))));
        Box box2;
        pp.get("b2", box2);
        AMREX_ALWAYS_ASSERT(box == box2);

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
    { // boolean strings queried as int
        ParmParse pp("bool");
        int v = -1;
        pp.get("true_val", v);
        AMREX_ALWAYS_ASSERT(v == 1);
        pp.get("false_val", v);
        AMREX_ALWAYS_ASSERT(v == 0);
        pp.get("True_val", v);
        AMREX_ALWAYS_ASSERT(v == 1);
        pp.get("FALSE_val", v);
        AMREX_ALWAYS_ASSERT(v == 0);
        pp.get("t_val", v);
        AMREX_ALWAYS_ASSERT(v == 1);
        pp.get("f_val", v);
        AMREX_ALWAYS_ASSERT(v == 0);
        long lv = -1;
        pp.get("true_val", lv);
        AMREX_ALWAYS_ASSERT(lv == 1);
        pp.get("false_val", lv);
        AMREX_ALWAYS_ASSERT(lv == 0);
        long long llv = -1;
        pp.get("true_val", llv);
        AMREX_ALWAYS_ASSERT(llv == 1);
        pp.get("false_val", llv);
        AMREX_ALWAYS_ASSERT(llv == 0);
    }
    { // boolean strings queried as bool
        ParmParse pp("bool");
        bool v = false;
        pp.get("true_val", v);
        AMREX_ALWAYS_ASSERT(v == true);
        pp.get("false_val", v);
        AMREX_ALWAYS_ASSERT(v == false);
        pp.get("True_val", v);
        AMREX_ALWAYS_ASSERT(v == true);
        pp.get("FALSE_val", v);
        AMREX_ALWAYS_ASSERT(v == false);
        pp.get("t_val", v);
        AMREX_ALWAYS_ASSERT(v == true);
        pp.get("f_val", v);
        AMREX_ALWAYS_ASSERT(v == false);
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
    { // AMREX_SPACEDIM
        ParmParse pp("macro.spacedim");
        std::vector<int> n_cell;
        double t;
        int use_gpu;
        pp.getarr("n_cell", n_cell);
        pp.get("t", t);
        pp.get("use_gpu", use_gpu);
        AMREX_ALWAYS_ASSERT(n_cell.size() == AMREX_SPACEDIM);
#if (AMREX_SPACEDIM == 1)
        AMREX_ALWAYS_ASSERT(n_cell[0] == 256);
#elif (AMREX_SPACEDIM == 2)
        AMREX_ALWAYS_ASSERT(n_cell[0] == 128 && n_cell[1] == 128);
#else
        AMREX_ALWAYS_ASSERT(n_cell[0] == 64 && n_cell[1] == 64 && n_cell[2] == 64);
#endif
#if (AMREX_SPACEDIM >= 2)
        AMREX_ALWAYS_ASSERT(almostEqual(t,0.5));
#else
        AMREX_ALWAYS_ASSERT(almostEqual(t,1.5));
#endif
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(use_gpu == 1);
#else
        AMREX_ALWAYS_ASSERT(use_gpu == 0);
#endif
    }
    { // AMREX_USE_GPU
        ParmParse pp("macro.use_gpu");
        int foo, spacedim;
        std::string bar;
        pp.get("foo", foo);
        pp.get("bar", bar);
        pp.get("spacedim", spacedim);
#ifdef AMREX_USE_GPU
        AMREX_ALWAYS_ASSERT(foo == 64 && bar == "use_gpu" && spacedim == AMREX_SPACEDIM*10);
#else
        AMREX_ALWAYS_ASSERT(foo == 32 && bar == "use_cpu" && spacedim == AMREX_SPACEDIM);
#endif
    }
    { // add & addarr
        ParmParse pp;
        pp.add("bt", true);
        pp.add("bf", false);
        std::string s;
        pp.get("bt", s); // It is intentional to read bool as string
        AMREX_ALWAYS_ASSERT(s == "true");
        pp.get("bf", s); // It is intentional to read bool as string
        AMREX_ALWAYS_ASSERT(s == "false");

        pp.add("doubleone",double(1));
        pp.get("doubleone", s); // It is intentional to read double as string
        AMREX_ALWAYS_ASSERT(s == "1.0");
        int intone;
        pp.query("doubleone", intone); // We are allowed to read 1.0 as 1.
        AMREX_ALWAYS_ASSERT(intone == 1);

        pp.add("string_scalar", "An string with white spaces");
        pp.get("string_scalar", s);
        AMREX_ALWAYS_ASSERT(s == "An string with white spaces");

        std::vector<std::string> sv{"string a", " string b", " string c ", "string-d"};
        pp.addarr("string_vector", sv);
        for (int i = 0; i < int(sv.size()); ++i) {
            pp.get("string_vector", s, i);
            AMREX_ALWAYS_ASSERT(s == sv[i]);
        }
    }
    if (ParallelDescriptor::IOProcessor()) // print & addfile
    {
        {
            ParmParse pp;
            pp.add("string-for-testing-addfile", "string for testing addfile");
            pp.add("string-for-testing-addfile", "string for testing addfile");
            int n = pp.countname("string-for-testing-addfile");
            AMREX_ALWAYS_ASSERT(n==2);
        }
        std::ofstream ofs("my-inputs");
        ParmParse::prettyPrintTable(ofs);
        ofs.close();
        ParmParse::addfile("my-inputs");
        std::string s;
        ParmParse pp;
        pp.get("string-for-testing-addfile", s);
        int n = pp.countname("string-for-testing-addfile");
        AMREX_ALWAYS_ASSERT(n==3 && s == "string for testing addfile");
    }
    {
        amrex::Print() << "SUCCESS\n";
    }
    amrex::Finalize();
}
