#include <AMReX_OpenMP.H>
#include <AMReX.H>
#include <AMReX_ParmParse.H>
#include <AMReX_Print.H>

#if defined(__APPLE__)
#include <sys/types.h>
#include <sys/sysctl.h>
#endif

#if defined(_WIN32)
#include <windows.h>
#endif

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>


namespace amrex
{
    int
    numUniquePhysicalCores ()
    {
        int ncores;

#if defined(__APPLE__)
        size_t len = sizeof(ncores);
        // See hw.physicalcpu and hw.physicalcpu_max
        //   https://developer.apple.com/documentation/kernel/1387446-sysctlbyname/determining_system_capabilities/
        //   https://developer.apple.com/documentation/kernel/1387446-sysctlbyname
        if (sysctlbyname("hw.physicalcpu", &ncores, &len, NULL, 0) == -1) {
            if (system::verbose > 0) {
                amrex::Print() << "numUniquePhysicalCores(): Error receiving hw.physicalcpu! "
                               << "Defaulting to visible cores.\n";
            }
            ncores = int(std::thread::hardware_concurrency());
        }
#elif defined(__linux__)
        std::set<std::vector<int>> uniqueThreadSets;
        int cpuIndex = 0;

        while (true) {
            // for each logical CPU in cpuIndex from 0...N-1
            std::string path = "/sys/devices/system/cpu/cpu" + std::to_string(cpuIndex) + "/topology/thread_siblings_list";
            std::ifstream file(path);
            if (!file.is_open()) {
                break; // no further CPUs to check
            }

            // find its siblings
            std::vector<int> siblings;
            std::string line;
            if (std::getline(file, line)) {
                std::stringstream ss(line);
                std::string token;

                // Possible syntax: 0-3, 8-11, 14,17
                // https://github.com/torvalds/linux/blob/v6.5/Documentation/ABI/stable/sysfs-devices-system-cpu#L68-L72
                while (std::getline(ss, token, ',')) {
                    size_t dashPos = token.find('-');
                    if (dashPos != std::string::npos) {
                        // Range detected
                        int start = std::stoi(token.substr(0, dashPos));
                        int end = std::stoi(token.substr(dashPos + 1));
                        for (int i = start; i <= end; ++i) {
                            siblings.push_back(i);
                        }
                    } else {
                        siblings.push_back(std::stoi(token));
                    }
                }
            }

            // and record the siblings group
            // (assumes: ascending and unique sets per cpuIndex)
            uniqueThreadSets.insert(siblings);
            cpuIndex++;
        }

        if (cpuIndex == 0) {
            if (system::verbose > 0) {
                amrex::Print() << "numUniquePhysicalCores(): Error reading CPU info.\n";
            }
            ncores = int(std::thread::hardware_concurrency());
        } else {
            ncores = int(uniqueThreadSets.size());
        }
#elif defined(_WIN32)
        DWORD length = 0;
        bool result = GetLogicalProcessorInformation(NULL, &length);

        if (!result) {
            if (system::verbose > 0) {
                amrex::Print() << "numUniquePhysicalCores(): Failed to get logical processor information! "
                               << "Defaulting to visible cores.\n";
            }
            ncores = int(std::thread::hardware_concurrency());
        }
        else {
            std::vector<SYSTEM_LOGICAL_PROCESSOR_INFORMATION> buffer(length / sizeof(SYSTEM_LOGICAL_PROCESSOR_INFORMATION));
            if (!GetLogicalProcessorInformation(&buffer[0], &length)) {
                if (system::verbose > 0) {
                    amrex::Print() << "numUniquePhysicalCores(): Failed to get logical processor information! "
                                   << "Defaulting to visible cores.\n";
                }
                ncores = int(std::thread::hardware_concurrency());
            } else {
                ncores = 0;
                for (const auto& info : buffer) {
                    if (info.Relationship == RelationProcessorCore) {
                        ncores++;
                    }
                }
            }
        }
#else
        // TODO:
        //   BSD
        if (system::verbose > 0) {
            amrex::Print() << "numUniquePhysicalCores(): Unknown system. Defaulting to visible cores.\n";
        }
        ncores = int(std::thread::hardware_concurrency());
#endif
        return ncores;
    }
} // namespace amrex

#ifdef AMREX_USE_OMP
namespace amrex::OpenMP
{
    namespace {
        constexpr int nlocks = 128;
        omp_lock_t omp_locks[nlocks];
        unsigned int initialized = 0;
    }

    void Initialize ()
    {
        if (initialized) {
            ++initialized;
            return;
        }

        amrex::ParmParse pp("amrex");
        std::string omp_threads = "system";
        pp.queryAdd("omp_threads", omp_threads);

        auto to_int = [](std::string const & str_omp_threads) {
            std::optional<int> num = std::stoi(str_omp_threads);
            return num;
        };

        if (omp_threads == "system") {
            // default or OMP_NUM_THREADS environment variable
        } else if (omp_threads == "nosmt") {
            char const *env_omp_num_threads = std::getenv("OMP_NUM_THREADS");
            if (env_omp_num_threads == nullptr) {
                omp_set_num_threads(numUniquePhysicalCores());
            }
            else if (amrex::system::verbose > 1) {
                amrex::Print() << "amrex.omp_threads was set to nosmt,"
                               << "but OMP_NUM_THREADS was set. Will keep "
                               << "OMP_NUM_THREADS=" << env_omp_num_threads << ".\n";
            }
        } else {
            std::optional<int> num_omp_threads = to_int(omp_threads);
            if (num_omp_threads.has_value()) {
                omp_set_num_threads(num_omp_threads.value());
            }
            else {
                if (amrex::system::verbose > 0) {
                    amrex::Print() << "amrex.omp_threads has an unknown value: "
                                   << omp_threads
                                   << " (try system, nosmt, or a positive integer)\n";
                }
            }
        }

        for (auto& lck : omp_locks) {
            omp_init_lock(&lck);
        }

        ++initialized;
    }

    void Finalize ()
    {
        if (initialized) {
            --initialized;
            if (initialized == 0) {
                for (auto& lck : omp_locks) {
                    omp_destroy_lock(&lck);
                }
            }
        }
    }

    omp_lock_t* get_lock (int ilock)
    {
        ilock = ilock % nlocks;
        if (ilock < 0) { ilock += nlocks; }
        return omp_locks + ilock;
    }

} // namespace amrex::OpenMP
#endif // AMREX_USE_OMP
