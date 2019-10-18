#include "QuantumSyncEngineTableBuilder.H"

//Include the full QS engine with table generation support
//(after some consistency tests)
#ifdef PXRMP_CORE_ONLY
    #error The Table Builder is incompatible with PXRMP_CORE_ONLY
#endif

#ifdef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
    #warning quantum_sync_engine.hpp should not have been included before reaching this point.
#endif
#include "quantum_sync_engine.hpp"
//_______________________________________________
