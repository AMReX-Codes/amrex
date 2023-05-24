#include <AMReX.H>
#include <AMReX_Algorithm.H>
#include <AMReX_Parser.H>

#include <regex>

using namespace amrex;

double f (int icase, double x, double y, double z);

bool test_parser(int icase, std::string const& expr)
{
    double x = 1.23, y = 2.34, z = 3.45;
    auto result_native = f(icase, x, y, z);

    Parser parser(expr);
    parser.registerVariables({"x","y","z"});
    auto const exe = parser.compile<3>();
    auto result_parser = exe(x,y,z);

    amrex::Print() << "\ncase " << icase << ": " << expr << "\n";
    parser.printExe();

    return amrex::almostEqual(result_native, result_parser, 10);
}

int main (int argc, char* argv[])
{
    amrex::Initialize(argc, argv, false);

    {
        std::ifstream ifs("fn.cpp");
        if (!ifs.is_open()) {
            amrex::Abort("Failed to open fn.cpp");
        }

        std::regex case_re("case[[:space:]]+([[:digit:]]+)[[:space:]]*:");
        std::regex expr_re("return[[:space:]]*(.+)[[:space:]]*;");

        int ntests = 0;
        Vector<std::pair<int,std::string>> failed_tests;

        std::string line;
        std::smatch cm, em;
        while (std::getline(ifs, line)) {
            if (std::regex_search(line, cm, case_re)) {
                int icase = std::stoi(cm[1]);
                if (std::getline(ifs, line)) {
                    std::regex_search(line, em, expr_re);
                    std::string expr(em[1]);
                    ++ntests;
                    if (! test_parser(icase, expr)) {
                        failed_tests.push_back(std::make_pair(icase,expr));
                    }
                } else {
                    amrex::Abort("How did this happend? No line after case.");
                }
            }
        }

        amrex::Print() << "\n";
        if (failed_tests.empty()) {
            amrex::Print() << "All " << ntests << " tests passed.\n";
        } else {
            amrex::Print() << failed_tests.size() << " out of " << ntests
                           << " tests failed.\n";
            for (auto const& ie : failed_tests) {
                amrex::Print() << "  case " << ie.first << ": " << ie.second << "\n";
            }
        }
    }

    amrex::Finalize();
}
