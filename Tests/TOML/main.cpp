#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>
#include <AMReX_String.H>

using namespace amrex;

int main(int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    // Let's add `amrex.command_foo=6` to argv
    std::vector<char*> argv_v;
    argv_v.reserve(argc+1);
    for (int i = 0; i < argc; ++i) {
        argv_v.push_back(strdup(argv[i]));
    }
    argv_v.push_back(strdup("amrex.command_foo=6"));

    int argc_2 = argc + 1;
    char** argv_2 = argv_v.data();
    amrex::Initialize(argc_2, argv_2);

    if (argc == 1) {
        ParmParse::addfile("config.toml");
    }

#if (AMREX_SPACEDIM != 2)
    // string
    {
        ParmParse ppa("a");
        std::string str, str2, str3, str4;
        ppa.get("str", str);
        ppa.get("str2", str2);
        ppa.get("str3", str3);
        ppa.get("str4", str4);
        std::string::size_type pos = str.find('\n', 0);
        if (pos != std::string::npos) {
            str.replace(pos, 1, "\\n");
        }
        AMREX_ALWAYS_ASSERT(str == str2);
        AMREX_ALWAYS_ASSERT(str3 == "Here is one quotation mark: \". Simple enough.");
        AMREX_ALWAYS_ASSERT(str4 == "Here are two quotation marks: \"\". Simple enough.");
    }

    // array
    {
        ParmParse pp;
        std::vector<int> iv;
        pp.getarr("integers", iv);
        int n = pp.countval("integers");
        AMREX_ALWAYS_ASSERT(n == 3 && iv == std::vector<int>({1,2,3}));

        std::vector<std::string> colors;
        pp.getarr("colors", colors);
        n = pp.countval("colors");
        AMREX_ALWAYS_ASSERT(n == 3 && colors == std::vector<std::string>({"red","yellow","green"}));

        std::vector<std::string> aos;
        pp.getarr("array_of_strings", aos);
        n = pp.countval("array_of_strings");
        AMREX_ALWAYS_ASSERT(n == 3 && aos == std::vector<std::string>({"face[0]","face[1]","face[2]"}));

        std::vector<std::string> aos2;
        pp.getarr("array_of_strings_2", aos2);
        n = pp.countval("array_of_strings_2");
        AMREX_ALWAYS_ASSERT(n == 4 && aos2 == std::vector<std::string>({"face[-1]","f,ace[0]","fa,ce[1]","face[2],"}));

        std::vector<std::vector<int>> ivv;
        pp.getarr("nested_arrays_of_ints", ivv);
        n = pp.countval("nested_arrays_of_ints"); // size of the outer array
        AMREX_ALWAYS_ASSERT(n == 2 && ivv == std::vector<std::vector<int>>({{1,2},{3,4,5}}));

        std::vector<std::vector<std::string>> svv;
        pp.getarr("nested_arrays_of_strings", svv);
        n = pp.countval("nested_arrays_of_strings");
        AMREX_ALWAYS_ASSERT(n == 2 && svv == std::vector<std::vector<std::string>>(
                                {{"aaa","bbb"},{"ccc"}}));

        std::vector<std::vector<std::string>> svv2;
        pp.getarr("nested_arrays_of_strings_2", svv2);
        n = pp.countval("nested_arrays_of_strings_2");
        AMREX_ALWAYS_ASSERT(n == 2 && svv2 == std::vector<std::vector<std::string>>(
                                {{"aaa","bbb"},{"ccc]"}}));

        std::vector<std::vector<std::string>> svv3;
        pp.getarr("nested_arrays_of_strings_3", svv3);
        n = pp.countval("nested_arrays_of_strings_3");
        AMREX_ALWAYS_ASSERT(n == 2 && svv3 == std::vector<std::vector<std::string>>(
                                {{"aa[a"," bbb"},{"c\\\"c\\\" c "}}));

        std::vector<std::vector<std::string>> svv4;
        pp.getarr("nested_arrays_of_strings_4", svv4);
        n = pp.countval("nested_arrays_of_strings_4");
        AMREX_ALWAYS_ASSERT(n == 2 && svv4 == std::vector<std::vector<std::string>>(
                                {{"aa]a"," b,bb"},{"ccc["}}));

        std::vector<int> iv2;
        pp.getarr("integers2", iv2);
        n = pp.countval("integers2");
        AMREX_ALWAYS_ASSERT(n == 3 && iv == iv2);

        std::vector<int> iv3;
        pp.getarr("integers3", iv3);
        n = pp.countval("integers3");
        AMREX_ALWAYS_ASSERT(n == 2 && iv3 == std::vector<int>({1,2}));

        // Test indexed access to TOML arrays via get("name", val, ival)
        int ival;
        pp.get("integers", ival, 0);
        AMREX_ALWAYS_ASSERT(ival == 1);
        pp.get("integers", ival, 1);
        AMREX_ALWAYS_ASSERT(ival == 2);
        pp.get("integers", ival, 2);
        AMREX_ALWAYS_ASSERT(ival == 3);

        // Test query variant
        int qval = -1;
        bool found = pp.query("integers", qval, 1);
        AMREX_ALWAYS_ASSERT(found && qval == 2);
    }
#endif

    // table
    {
        std::string key1;
        int key2;
        ParmParse pp;
        int foo_val = -1;
        pp.get("foo", foo_val);
        AMREX_ALWAYS_ASSERT(foo_val == 5);
        pp.get("table-1.key1", key1);
        pp.get("table-1.key2", key2);
        AMREX_ALWAYS_ASSERT(key1 == "some string" && key2 == 123);

        ParmParse pp2("table-2");
        pp2.get("key1", key1);
        pp2.get("key2", key2);
        AMREX_ALWAYS_ASSERT(key1 == "another string" && key2 == 456);

        int v;
        ParmParse ppxy("x.y");
        ppxy.get("z.w.v", v);
        AMREX_ALWAYS_ASSERT(v == 1);

        pp.get("x.v", v);
        AMREX_ALWAYS_ASSERT(v == 2);

        std::string color;
        bool sweet=false;
        bool smooth=false;
        ParmParse ppfruit("fruit");
        ppfruit.get("apple.color", color);
        ppfruit.get("apple.taste.sweet", sweet);
        ppfruit.get("apple.texture.smooth", smooth);
        AMREX_ALWAYS_ASSERT(color == "red" && sweet && smooth);
    }

    // FILE
    {
        ParmParse pp("config-2");
        int x, y, z;
        pp.get("x", x);
        pp.get("y", y);
        pp.get("z", z);
        AMREX_ALWAYS_ASSERT(x == 10 && y == 20 && z == 30);
    }

    // Error handling
    {
        ParmParse pp("invalid");
        std::vector<int> val;
        try {
            pp.getarr("empty_array", val);
            amrex::Abort("Should not get here, because getarr should raise an exception on an empty array");
        } catch (std::runtime_error const& e) {
            amrex::ignore_unused(e);
        }
        try {
            pp.getarr("invalid_array", val);
            amrex::Abort("Should not get here, because getarr should raise an exception on an invalid array");
        } catch (std::runtime_error const& e) {
            amrex::ignore_unused(e);
        }
    }

    // command line
    {
        ParmParse pp;
        int command_foo;
        pp.get("amrex.command_foo", command_foo);
        AMREX_ALWAYS_ASSERT(command_foo == 6);
    }

    // explicitly added
    {
        ParmParse pp("amrex");
        pp.add("bar", 137);
        int c;
        pp.get("bar", c);
        AMREX_ALWAYS_ASSERT(c == 137);
    }

    amrex::Finalize();

    for (auto* p : argv_v) {
        std::free(p);
    }

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
}
