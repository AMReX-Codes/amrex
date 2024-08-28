#ifndef JOKE_WARPX
#define JOKE_WARPX

namespace WarpX
{
  static short electromagnetic_solver_id = 2;
  static bool do_dive_cleaning = true;
  static bool use_filter = true;
  static amrex::IntVect filter_npass_each_dir = {11,22,33};
  static int n_rz_azimuthal_modes = 999;
  static int nox = 77;
  static int noy = 88;
  static int noz = 99;
  static short particle_pusher_algo = 0;
  static short field_gathering_algo = 1;
  static short current_deposition_algo = 0;

  struct CurrentDepositionAlgo {
    enum {
         Esirkepov = 0,
         Direct = 1,
         Vay = 2
    };
  };

  struct GatheringAlgo {
    enum {
         EnergyConserving = 0,
         MomentumConserving
    };
  };

  struct ParticlePusherAlgo {
    enum {
        Boris = 0,
        Vay = 1,
        HigueraCary = 2
    };
  };

  struct ElectromagneticSolverAlgo {
    enum {
      None = 0,
      Yee = 1,
      CKC = 2,
      PSATD = 3,
      Other = 4
    };
  };

}


#endif
