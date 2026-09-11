#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_ParmParse.H>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr char warning[] =
    "TinyProfiler warning: values marked with ? are asynchronous GPU timings and may be inaccurate.";

void generate_profile ()
{
    {
        BL_PROFILE("ordinary_only");
    }
    {
        BL_PROFILE_ASYNC("async_scoped");
    }

    BL_PROFILE_VAR_ASYNC("async_variable", async_variable);
    BL_PROFILE_VAR_STOP(async_variable);

    BL_PROFILE_VAR_NS_ASYNC("async_no_start", async_no_start);
    BL_PROFILE_VAR_START(async_no_start);
    BL_PROFILE_VAR_STOP(async_no_start);

    {
        BL_PROFILE("same_name");
    }
    {
        BL_PROFILE_ASYNC("same_name");
    }

    if (amrex::ParallelDescriptor::NProcs() == 1 ||
        amrex::ParallelDescriptor::MyProc() == 1) {
        BL_PROFILE_ASYNC("rank_local_async");
    }

    {
        BL_PROFILE_REGION("ordinary_region");
        {
            BL_PROFILE("ordinary_in_ordinary_region");
        }
        {
            BL_PROFILE_ASYNC("async_in_ordinary_region");
        }
    }

    {
        BL_PROFILE_REGION_ASYNC("async_region");
        {
            BL_PROFILE("ordinary_in_async_region");
        }
        {
            BL_PROFILE_ASYNC("async_in_async_region");
        }
    }

    {
        BL_PROFILE_REGION("same_region");
        BL_PROFILE("ordinary_same_region_timer");
    }
    {
        BL_PROFILE_REGION_ASYNC("same_region");
        BL_PROFILE("async_same_region_timer");
    }

    if (amrex::ParallelDescriptor::NProcs() == 1 ||
        amrex::ParallelDescriptor::MyProc() == 1) {
        BL_PROFILE_REGION_ASYNC("rank_local_async_region");
        BL_PROFILE("ordinary_rank_local_region_timer");
    }
}

bool contains (std::string const& text, std::string const& value)
{
    return text.find(value) != std::string::npos;
}

int count (std::string const& text, std::string const& value)
{
    int result = 0;
    std::size_t pos = 0;
    while ((pos = text.find(value, pos)) != std::string::npos) {
        ++result;
        pos += value.size();
    }
    return result;
}

std::vector<std::string> rows_named (std::string const& table, std::string const& name)
{
    std::vector<std::string> result;
    std::istringstream lines(table);
    for (std::string line; std::getline(lines, line); ) {
        std::istringstream fields(line);
        std::string row_name;
        if ((fields >> row_name) && row_name == name) {
            result.push_back(std::move(line));
        }
    }
    return result;
}

std::vector<std::string> data_rows (std::string const& table)
{
    std::vector<std::string> result;
    std::istringstream lines(table);
    for (std::string line; std::getline(lines, line); ) {
        std::istringstream fields(line);
        std::string name;
        std::string calls;
        std::string tmin;
        std::string tavg;
        std::string tmax;
        std::string percent;
        if ((fields >> name >> calls >> tmin >> tavg >> tmax >> percent) &&
            name != "Name") {
            result.push_back(std::move(line));
        }
    }
    return result;
}

std::vector<std::string> region_tables (std::string const& report, std::string const& name)
{
    std::vector<std::string> result;
    std::string const begin_marker = "BEGIN REGION " + name + "\n";
    std::string const end_marker = "END REGION " + name + "\n";
    std::size_t begin = 0;
    while ((begin = report.find(begin_marker, begin)) != std::string::npos) {
        std::size_t const end = report.find(end_marker, begin);
        if (end == std::string::npos) {
            break;
        }
        result.push_back(report.substr(begin, end + end_marker.size() - begin));
        begin = end + end_marker.size();
    }
    return result;
}

bool marked_row (std::string const& row)
{
    std::istringstream fields(row);
    std::string name;
    std::string calls;
    std::string tmin;
    std::string tavg;
    std::string tmax;
    std::string percent;
    if (!(fields >> name >> calls >> tmin >> tavg >> tmax >> percent)) {
        return false;
    }
    auto ends_with = [] (std::string const& value, std::string const& suffix) {
        return value.size() >= suffix.size() &&
               value.compare(value.size()-suffix.size(), suffix.size(), suffix) == 0;
    };
    return !contains(name, "?") && !contains(calls, "?") &&
           ends_with(tmin, "?") && ends_with(tavg, "?") &&
           ends_with(tmax, "?") && ends_with(percent, "?%");
}

bool unmarked_row (std::string const& row)
{
    return !contains(row, "?");
}

long calls_in_row (std::string const& row)
{
    std::istringstream fields(row);
    std::string name;
    long calls = -1;
    fields >> name >> calls;
    return calls;
}

bool expect (bool condition, std::string const& message)
{
    if (!condition) {
        std::cerr << "TinyProfiler async test failure: " << message << '\n';
    }
    return condition;
}

bool validate_tiny_report (std::string const& output, bool hide_async, bool expect_other)
{
    bool ok = true;
    std::size_t const first_region = output.find("BEGIN REGION");
    std::string const main_table = output.substr(0, first_region);

#ifdef AMREX_USE_GPU
    if (expect_other) {
        auto const other = rows_named(main_table, "Other");
        ok = expect(other.size() == 2 && marked_row(other[0]) && marked_row(other[1]),
                    "Other did not retain its async classification") && ok;
        ok = expect(count(output, warning) == 1,
                    "Other report must contain exactly one async warning") && ok;
        return ok;
    }

    if (hide_async) {
        ok = expect(!contains(output, "?"), "hidden report contains an async marker") && ok;
        ok = expect(!contains(output, warning), "hidden report contains the async warning") && ok;
        ok = expect(!contains(main_table, "async_scoped") &&
                    !contains(main_table, "async_variable") &&
                    !contains(main_table, "async_no_start") &&
                    !contains(main_table, "rank_local_async"),
                    "hidden report retained an async timer") && ok;
        ok = expect(region_tables(output, "async_region").empty() &&
                    region_tables(output, "rank_local_async_region").empty(),
                    "hidden report retained an async region table") && ok;

        auto const same_name = rows_named(main_table, "same_name");
        ok = expect(same_name.size() == 2 && calls_in_row(same_name[0]) == 1 &&
                    calls_in_row(same_name[1]) == 1,
                    "hidden report did not preserve the ordinary same-name timer") && ok;
        ok = expect(rows_named(main_table, "ordinary_in_async_region").size() == 2,
                    "hidden report removed an ordinary timer nested in an async region") && ok;

        auto const ordinary_regions = region_tables(output, "ordinary_region");
        ok = expect(ordinary_regions.size() == 1 &&
                    rows_named(ordinary_regions[0], "async_in_ordinary_region").empty() &&
                    rows_named(ordinary_regions[0], "ordinary_in_ordinary_region").size() == 2,
                    "hidden report filtered an ordinary region incorrectly") && ok;
        ok = expect(region_tables(output, "same_region").size() == 1,
                    "hidden report did not preserve only the ordinary same-name region") && ok;
    } else {
        ok = expect(count(output, warning) == 1,
                    "default report must contain exactly one async warning") && ok;

        for (auto const& name : {"async_scoped", "async_variable", "async_no_start",
                                 "rank_local_async", "async_in_ordinary_region"}) {
            auto const rows = rows_named(main_table, name);
            ok = expect(rows.size() == 2 && marked_row(rows[0]) && marked_row(rows[1]),
                        std::string("default report did not mark ") + name) && ok;
        }

        auto const same_name = rows_named(main_table, "same_name");
        int marked = 0;
        int unmarked = 0;
        for (auto const& row : same_name) {
            marked += marked_row(row) ? 1 : 0;
            unmarked += unmarked_row(row) ? 1 : 0;
        }
        ok = expect(same_name.size() == 4 && marked == 2 && unmarked == 2,
                    "ordinary and async same-name timers were not kept separate") && ok;

        auto const ordinary_in_async = rows_named(main_table, "ordinary_in_async_region");
        ok = expect(ordinary_in_async.size() == 2 &&
                    unmarked_row(ordinary_in_async[0]) && unmarked_row(ordinary_in_async[1]),
                    "async region changed a nested ordinary timer in the main table") && ok;

        auto const async_regions = region_tables(output, "async_region");
        ok = expect(async_regions.size() == 1,
                    "default report omitted the async region") && ok;
        if (async_regions.size() == 1) {
            auto const rows = data_rows(async_regions[0]);
            bool all_marked = !rows.empty();
            for (auto const& row : rows) {
                all_marked = all_marked && marked_row(row);
            }
            ok = expect(all_marked, "async region did not mark every timing row") && ok;
        }

        auto const same_regions = region_tables(output, "same_region");
        int marked_tables = 0;
        int unmarked_tables = 0;
        for (auto const& table : same_regions) {
            auto const rows = data_rows(table);
            bool all_marked = !rows.empty();
            bool all_unmarked = !rows.empty();
            for (auto const& row : rows) {
                all_marked = all_marked && marked_row(row);
                all_unmarked = all_unmarked && unmarked_row(row);
            }
            marked_tables += all_marked ? 1 : 0;
            unmarked_tables += all_unmarked ? 1 : 0;
        }
        ok = expect(same_regions.size() == 2 && marked_tables == 1 && unmarked_tables == 1,
                    "ordinary and async same-name regions were not kept separate") && ok;
    }
#else
    ok = expect(!contains(output, "?") && !contains(output, warning),
                "CPU aliases produced async annotations") && ok;
    if (expect_other) {
        auto const other = rows_named(main_table, "Other");
        ok = expect(other.size() == 2 && unmarked_row(other[0]) && unmarked_row(other[1]),
                    "CPU Other row did not remain ordinary") && ok;
        return ok;
    }

    ok = expect(contains(main_table, "async_scoped") &&
                contains(main_table, "async_variable") &&
                contains(main_table, "async_no_start") &&
                contains(main_table, "rank_local_async"),
                "CPU aliases did not produce ordinary timer output") && ok;

    auto const same_name = rows_named(main_table, "same_name");
    ok = expect(same_name.size() == 2 && calls_in_row(same_name[0]) == 2 &&
                calls_in_row(same_name[1]) == 2,
                "CPU async alias did not aggregate with its ordinary same-name timer") && ok;
    ok = expect(region_tables(output, "async_region").size() == 1 &&
                region_tables(output, "same_region").size() == 1,
                "CPU async region alias did not behave as an ordinary region") && ok;
#endif

    return ok;
}

}

int main (int argc, char* argv[])
{
#ifdef AMREX_USE_MPI
    MPI_Init(&argc, &argv);
#endif

    amrex::Initialize(argc, argv);

    bool hide_async = false;
    bool expect_other = false;
    std::string output_file;
    amrex::ParmParse("test").query("hide_async", hide_async);
    amrex::ParmParse("test").query("expect_other", expect_other);
    amrex::ParmParse("tiny_profiler").query("output_file", output_file);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(output_file.c_str());
    }
    amrex::ParallelDescriptor::Barrier();

    generate_profile();

    bool const io_processor = amrex::ParallelDescriptor::IOProcessor();
    amrex::Finalize();

    int result = 0;
#ifdef AMREX_TINY_PROFILING
    if (io_processor) {
        std::ifstream stream(output_file);
        std::string const output((std::istreambuf_iterator<char>(stream)),
                                 std::istreambuf_iterator<char>());
        if (!expect(stream.good() || stream.eof(), "could not read profiler output") ||
            !validate_tiny_report(output, hide_async, expect_other)) {
            result = 1;
        }
    }
#else
    amrex::ignore_unused(io_processor, hide_async, expect_other);
#endif

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
