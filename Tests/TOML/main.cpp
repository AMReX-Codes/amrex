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
        AMREX_ALWAYS_ASSERT(iv == std::vector<int>({1,2,3}));

        std::vector<std::string> colors;
        pp.getarr("colors", colors);
        AMREX_ALWAYS_ASSERT(colors == std::vector<std::string>({"red","yellow","green"}));

        std::vector<std::vector<int>> ivv;
        pp.getarr("nested_arrays_of_ints", ivv);
        AMREX_ALWAYS_ASSERT(ivv == std::vector<std::vector<int>>({{1,2},{3,4,5}}));

        std::vector<std::vector<std::string>> svv;
        pp.getarr("nested_arrays_of_strings", svv);
        AMREX_ALWAYS_ASSERT(svv == std::vector<std::vector<std::string>>(
                                {{"aaa","bbb"},{"ccc"}}));

        std::vector<int> iv2;
        pp.getarr("integers2", iv2);
        AMREX_ALWAYS_ASSERT(iv == iv2);

        std::vector<int> iv3;
        pp.getarr("integers3", iv3);
        AMREX_ALWAYS_ASSERT(iv3 == std::vector<int>({1,2}));
    }
#endif

    // table
    {
        std::string key1;
        int key2;
        ParmParse pp;
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
