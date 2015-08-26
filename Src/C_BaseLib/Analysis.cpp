#include <Analysis.H>

namespace Analysis
{
    // Various signals to send to the sidecar group.
    const int NyxHaloFinderSignal = 42;
    const int QuitSignal = -1;

    AnalysisContainer *analysis;
}
