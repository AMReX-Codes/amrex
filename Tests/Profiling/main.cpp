#include <AMReX.H>
#include <AMReX_BLProfiler.H>
#include <AMReX_GpuContainers.H>

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr char output_file[] = "tiny_profiler_test.out";

void preactivation_phase ()
{
    BL_PROFILE("preactivation_outermost");

    {
        BL_PROFILE("preactivation_nested");
    }

    {
        BL_PROFILE_REGION("preactivation_region");
        BL_PROFILE("preactivation_region_timer");
    }
}

void collective_activation_phase ()
{
    BL_PROFILE("ordinary_outermost");

    {
        BL_PROFILE("ordinary_nested");
    }

    if (amrex::ParallelDescriptor::NProcs() == 1 ||
        amrex::ParallelDescriptor::MyProc() == 1) {
        BL_PROFILE_GPU_SYNC("rank_one_gpu_sync");
    }

    {
        BL_PROFILE_REGION("ordinary_region");
        BL_PROFILE("ordinary_region_timer");
    }
}

void macro_coverage_phase ()
{
    BL_PROFILE("coverage_outermost");

    {
        BL_PROFILE_GPU_SYNC("scoped_gpu_sync");
    }

    BL_PROFILE_VAR_GPU_SYNC("variable_gpu_sync", gpu_sync_var);
    BL_PROFILE_VAR_STOP(gpu_sync_var);

    BL_PROFILE_VAR_NS_GPU_SYNC("no_start_gpu_sync", no_start_gpu_sync_var);
    BL_PROFILE_VAR_START(no_start_gpu_sync_var);
    BL_PROFILE_VAR_STOP(no_start_gpu_sync_var);

    {
        BL_PROFILE("same_name");
    }
    {
        BL_PROFILE_GPU_SYNC("same_name");
    }

    {
        BL_PROFILE_REGION_GPU_SYNC("gpu_sync_region");
        {
            BL_PROFILE("ordinary_timer_in_gpu_sync_region");
        }
        {
            BL_PROFILE_GPU_SYNC("gpu_sync_timer_in_gpu_sync_region");
        }
    }

#if defined(AMREX_USE_GPU) && defined(AMREX_TINY_PROFILING)
    amrex::Gpu::DeviceVector<int> device_value(1);
    amrex::Gpu::PinnedVector<int> host_value(1, 0);
    int* device_ptr = device_value.data();

    {
        BL_PROFILE_GPU_SYNC("async_gpu_work");
        amrex::ParallelFor(1, [=] AMREX_GPU_DEVICE (int) noexcept {
            device_ptr[0] = 42;
        });
        amrex::Gpu::dtoh_memcpy_async(host_value.data(), device_ptr, sizeof(int));
    }

    AMREX_ALWAYS_ASSERT_WITH_MESSAGE(host_value[0] == 42,
        "BL_PROFILE_GPU_SYNC returned before work on the current GPU stream completed");
#endif
}

std::vector<std::string> split_reports (std::string const& output)
{
    constexpr char marker[] = "TinyProfiler total time across processes";
    std::vector<std::string> reports;
    std::size_t begin = output.find(marker);
    while (begin != std::string::npos) {
        std::size_t const next = output.find(marker, begin + 1);
        reports.emplace_back(output.substr(begin, next - begin));
        begin = next;
    }
    return reports;
}

bool contains (std::string const& text, std::string const& value)
{
    return text.find(value) != std::string::npos;
}

std::vector<long> calls_in_main_table (std::string const& report, std::string const& name)
{
    std::size_t const region_begin = report.find("BEGIN REGION");
    std::string const main_table = report.substr(0, region_begin);
    std::istringstream lines(main_table);
    std::vector<long> result;
    for (std::string line; std::getline(lines, line); ) {
        std::istringstream fields(line);
        std::string row_name;
        long calls = -1;
        if ((fields >> row_name >> calls) && row_name == name) {
            result.push_back(calls);
        }
    }
    return result;
}

bool expect (bool condition, std::string const& message)
{
    if (!condition) {
        std::cerr << "TinyProfiler test failure: " << message << '\n';
    }
    return condition;
}

bool validate_report ()
{
    std::ifstream stream(output_file);
    std::string const output((std::istreambuf_iterator<char>(stream)),
                             std::istreambuf_iterator<char>());
    auto const reports = split_reports(output);
    bool ok = expect(stream.good() || stream.eof(), "could not read profiler output") &&
              expect(reports.size() == 3, "expected two flush reports and one final report");
    if (reports.size() != 3) {
        return false;
    }

    auto const& preactivation_report = reports[0];
    auto const& activated_flush_report = reports[1];
    auto const& final_report = reports[2];

    ok = expect(contains(preactivation_report, "preactivation_nested") &&
                contains(preactivation_report, "BEGIN REGION preactivation_region"),
                "flush before activation did not preserve the complete report") && ok;

#ifdef AMREX_USE_GPU
    // Only rank 1 uses a synchronized macro before the flush.  Every rank must
    // nevertheless choose the focused report collectively.
    ok = expect(contains(activated_flush_report, "ordinary_outermost"),
                "focused flush omitted the outermost timer") && ok;
    ok = expect(contains(activated_flush_report, "rank_one_gpu_sync"),
                "focused flush omitted the rank-local synchronized timer") && ok;
    ok = expect(!contains(activated_flush_report, "ordinary_nested"),
                "focused flush retained an ordinary nested timer") && ok;
    ok = expect(!contains(activated_flush_report, "BEGIN REGION ordinary_region"),
                "focused flush retained an ordinary region") && ok;

    ok = expect(contains(final_report, "coverage_outermost"),
                "focused report omitted an outermost timer") && ok;
    ok = expect(contains(final_report, "scoped_gpu_sync") &&
                contains(final_report, "variable_gpu_sync") &&
                contains(final_report, "no_start_gpu_sync"),
                "focused report omitted a synchronized timer macro") && ok;
    ok = expect(contains(final_report, "BEGIN REGION gpu_sync_region"),
                "focused report omitted the synchronized region") && ok;
    ok = expect(contains(final_report, "gpu_sync_timer_in_gpu_sync_region"),
                "synchronized region omitted its synchronized timer") && ok;
    ok = expect(!contains(final_report, "ordinary_timer_in_gpu_sync_region"),
                "synchronized region retained an ordinary nested timer") && ok;

    auto const same_name_calls = calls_in_main_table(final_report, "same_name");
    ok = expect(same_name_calls.size() == 2 && same_name_calls[0] == 1 &&
                same_name_calls[1] == 1,
                "ordinary samples contaminated the identically named synchronized timer") && ok;
#else
    // On CPU, every new macro is an exact alias for its ordinary counterpart.
    ok = expect(contains(activated_flush_report, "ordinary_nested") &&
                contains(activated_flush_report, "BEGIN REGION ordinary_region"),
                "CPU synchronized aliases changed ordinary TinyProfiler output") && ok;
    ok = expect(contains(final_report, "ordinary_timer_in_gpu_sync_region") &&
                contains(final_report, "BEGIN REGION gpu_sync_region"),
                "CPU synchronized region alias changed ordinary region output") && ok;

    auto const same_name_calls = calls_in_main_table(final_report, "same_name");
    ok = expect(same_name_calls.size() == 2 && same_name_calls[0] == 2 &&
                same_name_calls[1] == 2,
                "CPU synchronized timer alias did not aggregate like an ordinary timer") && ok;
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

    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::remove(output_file);
    }
    amrex::ParallelDescriptor::Barrier();

    preactivation_phase();
    BL_PROFILE_TINY_FLUSH();
    collective_activation_phase();
    BL_PROFILE_TINY_FLUSH();
    macro_coverage_phase();

    bool const io_processor = amrex::ParallelDescriptor::IOProcessor();
    amrex::Finalize();

    int result = 0;
#ifdef AMREX_TINY_PROFILING
    if (io_processor && !validate_report()) {
        result = 1;
    }
#endif

#ifdef AMREX_USE_MPI
    MPI_Finalize();
#endif
    return result;
}
