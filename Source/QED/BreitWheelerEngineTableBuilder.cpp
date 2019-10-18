#include "BreitWheelerEngineTableBuilder.H"


//Include the full BW engine with table generation support
//(after some consistency tests)
#ifdef PXRMP_CORE_ONLY
    #error The Table Builder is incompatible with PXRMP_CORE_ONLY
#endif

#ifdef __PICSAR_MULTIPHYSICS_BREIT_WHEELER_ENGINE__
    #warning breit_wheeler_engine.hpp should not have been included before reaching this point.
#endif
#include "breit_wheeler_engine.hpp"
//_______________________________________________
