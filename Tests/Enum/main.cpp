#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

using namespace amrex;

AMREX_ENUM(MyColor,  red,  green, blue );

namespace my_namespace {
    AMREX_ENUM(MyColor, orange,     yellow,cyan );
}

AMREX_ENUM(MyColor2,
           red,       // 0
           chi=red,   // 0
           green,     // 1
           blue,      // 2
           Default = green); // 1

AMREX_ENUM(Location, entry = 1, exit = -1, after_exit);

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv);
    {
        auto const& names = amrex::getEnumNameStrings<MyColor>();
        auto const& names2 = amrex::getEnumNameStrings<my_namespace::MyColor>();
        amrex::Print() << "colors:";
        for (auto const& name : names) {
            amrex::Print() << " " << name;
        }
        amrex::Print() << "\n";
        amrex::Print() << "colors:";
        for (auto const& name : names2) {
            amrex::Print() << " " << name;
        }
        amrex::Print() << "\n";

        ParmParse pp;
        {
            auto color = static_cast<MyColor>(999);
            pp.query("color1", color);
            amrex::Print() << "color = " << amrex::getEnumNameString(color) << '\n';
            AMREX_ALWAYS_ASSERT(color == MyColor::red);
        }
        {
            auto color = static_cast<MyColor>(999);
            pp.get("color2", color);
            amrex::Print() << "color = " << amrex::getEnumNameString(color) << '\n';
            AMREX_ALWAYS_ASSERT(color == MyColor::green);
        }
        {
            auto color = static_cast<MyColor>(999);
            pp.get("color3", color);
            amrex::Print() << "color = " << amrex::getEnumNameString(color) << '\n';
            AMREX_ALWAYS_ASSERT(color == MyColor::blue);
        }
        {
            auto color = static_cast<MyColor>(999);
            try {
                pp.query("color4", color);
            } catch (std::runtime_error const& e) {
                amrex::Print() << "As expected, " << e.what() << '\n';
            }
            AMREX_ALWAYS_ASSERT(color == static_cast<MyColor>(999));
            try {
                pp.get_enum_case_insensitive("color4", color);
            } catch (std::runtime_error const& e) {
                amrex::Print() << "As expected, " << e.what() << '\n';
            }
            AMREX_ALWAYS_ASSERT(color == static_cast<MyColor>(999));
        }
        {
            auto color = static_cast<MyColor>(999);
            try {
                pp.query("color5", color);
            } catch (std::runtime_error const& e) {
                amrex::Print() << "As expected, " << e.what() << '\n';
            }
            AMREX_ALWAYS_ASSERT(color == static_cast<MyColor>(999));
            pp.query_enum_case_insensitive("color5", color);
            amrex::Print() << "color = " << amrex::getEnumNameString(color) << '\n';
            AMREX_ALWAYS_ASSERT(color == MyColor::blue);
        }
        {
            std::vector<my_namespace::MyColor> color;
            pp.getarr("colors", color);
            AMREX_ALWAYS_ASSERT(color.size() == 3 &&
                                color[0] == my_namespace::MyColor::cyan &&
                                color[1] == my_namespace::MyColor::yellow &&
                                color[2] == my_namespace::MyColor::orange);
            std::vector<my_namespace::MyColor> color2;
            pp.queryarr("colors", color2);
            AMREX_ALWAYS_ASSERT(color.size() == 3 &&
                                color == color2 &&
                                color[0] == my_namespace::MyColor::cyan &&
                                color[1] == my_namespace::MyColor::yellow &&
                                color[2] == my_namespace::MyColor::orange);
            amrex::Print() << "colors:";
            for (auto const& c : color) {
                amrex::Print() << " " << amrex::getEnumNameString(c);
            }
            amrex::Print() << "\n";
        }
    }
    {
        auto names2 = amrex::getEnumNameStrings<MyColor2>();
        amrex::Print() << "Names in " << amrex::getEnumClassName<MyColor2>() << "\n";
        for (auto const& name : names2) {
            amrex::Print() << "  " << name << "\n";
        }

        auto const& kv = amrex::getEnumNameValuePairs<MyColor2>();
        amrex::Print() << "Name : Value\n";
        for (auto const& item : kv) {
            amrex::Print() << "  " << item.first << ": "
                           << static_cast<int>(item.second) << "\n";
        }

        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(MyColor2::red) == "red");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(MyColor2::chi) == "red");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(MyColor2::green) == "green");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(MyColor2::blue) == "blue");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(MyColor2::Default) == "green");

        static_assert(MyColor2::red == MyColor2::chi);
        static_assert(MyColor2::green == MyColor2::Default);

        AMREX_ALWAYS_ASSERT(amrex::getEnum<MyColor2>("red") == MyColor2::red);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<MyColor2>("chi") == MyColor2::chi);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<MyColor2>("green") == MyColor2::green);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<MyColor2>("blue") == MyColor2::blue);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<MyColor2>("Default") == MyColor2::Default);

        auto my_blue = MyColor2::Default;
        ParmParse pp;
        pp.query_enum_case_insensitive("color5", my_blue);
        AMREX_ALWAYS_ASSERT(my_blue == MyColor2::blue);

        auto my_green = MyColor2::red;
        pp.query_enum_sloppy("color6", my_green, "-.");
        AMREX_ALWAYS_ASSERT(my_green == MyColor2::green);
    }


    {
        auto const& kv = amrex::getEnumNameValuePairs<Location>();
        amrex::Print() << "Name : Value\n";
        for (auto const& item : kv) {
            amrex::Print() << "  " << item.first << ": "
                           << static_cast<int>(item.second) << "\n";
        }
        AMREX_ALWAYS_ASSERT(static_cast<int>(Location::entry) == 1);
        AMREX_ALWAYS_ASSERT(static_cast<int>(Location::exit) == -1);
        AMREX_ALWAYS_ASSERT(static_cast<int>(Location::after_exit) == 0);
        AMREX_ALWAYS_ASSERT(amrex::toUnderlying(Location::entry) == 1);
        AMREX_ALWAYS_ASSERT(amrex::toUnderlying(Location::exit) == -1);
        AMREX_ALWAYS_ASSERT(amrex::toUnderlying(Location::after_exit) == 0);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<Location>("entry") == Location::entry);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<Location>("exit") == Location::exit);
        AMREX_ALWAYS_ASSERT(amrex::getEnum<Location>("after_exit") == Location::after_exit);
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(Location::entry) == "entry");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(Location::exit) == "exit");
        AMREX_ALWAYS_ASSERT(amrex::getEnumNameString(Location::after_exit) == "after_exit");
    }

    amrex::Finalize();
}
