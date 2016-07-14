#include <bl_random_c.H>

#include <set>
#include <limits>
#include <cstdint>

namespace 
{
    std::uint_fast32_t bl_rng_parallel_seed (int s, int rank, int nprocs)
    {
	if (rank < 0 || rank >= nprocs) {
	    std::cerr << "bl_rng_parallel_seed: " << rank << ", " << nprocs << std::endl;
	    backtrace_handler(6);
	}

	std::mt19937 eng(s);
	eng.discard(10000);

	std::uniform_int_distribution<std::uint_fast32_t> dist
	    (std::numeric_limits<std::uint_fast32_t>::min(),
	     std::numeric_limits<std::uint_fast32_t>::max());

	std::uint_fast32_t r;
	std::set<std::uint_fast32_t> seeds;

	while (seeds.size() != rank+1) {
	    r = dist(eng);
	    seeds.insert(r);
	};

	return r;
    }
}

extern "C"
{
    int bl_rng_random_uint_c ()
    {
	std::random_device rd;
	std::uniform_int_distribution<int> dist(0, std::numeric_limits<int>::max());
	return dist(rd);
    }

    // 
    // engine
    //
    using BLRngEngine = std::mt19937;

    void bl_rng_new_engine_c (BLRngEngine*& eng, int s, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	eng = new BLRngEngine(seed);
	eng->discard(1000000);  // warm up
    }

    void bl_rng_delete_engine_c (BLRngEngine* eng)
    {
	delete eng;
    }

    void bl_rng_save_engine_c (const BLRngEngine* eng, const char* name)
    {
	BLRng_save(*eng, name);
    }

    void bl_rng_restore_engine_c (BLRngEngine*& eng, const char* name)
    {
	eng = new BLRngEngine();
	BLRng_restore(*eng, name);
    }

    //
    // uniform real distribution
    //
    using BLRngUniformReal = std::uniform_real_distribution<double>;

    void bl_rng_new_uniform_real_c (BLRngUniformReal*& distro, double a, double b)
    {
	distro = new BLRngUniformReal(a,b);
    }

    void bl_rng_delete_uniform_real_c (BLRngUniformReal* distro)
    {
	delete distro;
    }

    double bl_rng_get_uniform_real_c (BLRngUniformReal* distro, BLRngEngine* eng)
    {
	return BLRng_get<double>(*distro,*eng);
    }

    void bl_rng_save_uniform_real_c (const BLRngUniformReal* distro, const char* name)
    {
	BLRng_save(*distro, name);
    }

    void bl_rng_restore_uniform_real_c (BLRngUniformReal*& distro, const char* name)
    {
	distro = new BLRngUniformReal();
	BLRng_restore(*distro, name);
    }

    //
    // normal distribution
    //
    using BLRngNormal = std::normal_distribution<double>;

    void bl_rng_new_normal_c (BLRngNormal*& distro, double mean, double stddev)
    {
	distro = new BLRngNormal(mean,stddev);
    }

    void bl_rng_delete_normal_c (BLRngNormal* distro)
    {
	delete distro;
    }

    double bl_rng_get_normal_c (BLRngNormal* distro, BLRngEngine* eng)
    {
	return BLRng_get<double>(*distro, *eng);
    }

    void bl_rng_save_normal_c (const BLRngNormal* distro, const char* name)
    {
	BLRng_save(*distro, name);
    }

    void bl_rng_restore_normal_c (BLRngNormal*& distro, const char* name)
    {
	distro = new BLRngNormal();
	BLRng_restore(*distro, name);
    }

    //
    // poisson distribution
    //
    using BLRngPoisson = std::poisson_distribution<int>;

    void bl_rng_new_poisson_c (BLRngPoisson*& distro, double mean)
    {
	distro = new BLRngPoisson(mean);
    }

    void bl_rng_delete_poisson_c (BLRngPoisson* distro)
    {
	delete distro;
    }

    int bl_rng_get_poisson_c (BLRngPoisson* distro, BLRngEngine* eng)
    {
	return BLRng_get<int>(*distro, *eng);
    }

    void bl_rng_save_poisson_c (const BLRngPoisson* distro, const char* name)
    {
	BLRng_save(*distro, name);
    }

    void bl_rng_restore_poisson_c (BLRngPoisson*& distro, const char* name)
    {
	distro = new BLRngPoisson();
	BLRng_restore(*distro, name);
    }

    //
    // binomial distribution
    //
    using BLRngBinomial = std::binomial_distribution<int>;

    void bl_rng_new_binomial_c (BLRngBinomial*& distro, int t, double p)
    {
	distro = new BLRngBinomial(t,p);
    }

    void bl_rng_delete_binomial_c (BLRngBinomial* distro)
    {
	delete distro;
    }

    int bl_rng_get_binomial_c (BLRngBinomial* distro, BLRngEngine* eng)
    {
	return BLRng_get<int>(*distro, *eng);
    }

    void bl_rng_save_binomial_c (const BLRngBinomial* distro, const char* name)
    {
	BLRng_save(*distro, name);
    }

    void bl_rng_restore_binomial_c (BLRngBinomial*& distro, const char* name)
    {
	distro = new BLRngBinomial();
	BLRng_restore(*distro, name);
    }

    void hg_genrand (double* rn, BLRngEngine* rng)
    {
	constexpr double fac = (1.0 - std::numeric_limits<double>::epsilon()) / std::mt19937::max();

	// compile time check!
	static_assert(std::mt19937::min() == 0, "hg_genrand: std::mt19937::min() != 0");
	static_assert((double)std::mt19937::max() * fac < 1.0, "hg_genrand: < 1 failed");
	// Some codes have the following constant hard wired.
	static_assert((double)UINT32_MAX * 2.3283064370807969e-10 < 1.0, 
		      "hg_genrand: check constant");
	static_assert((double)UINT32_MAX * 2.3283064370807969e-10 > 1.0 - 2.0*std::numeric_limits<double>::epsilon(),
		      "hg_genrand: check constant");

	auto y = (*rng)();
	*rn = (double)y * fac; /* reals: [0,1)-interval */    
    }

    void hg_genrand_sp (float* rn, BLRngEngine* rng)
    {    
	constexpr float fac = (1.0f - std::numeric_limits<float>::epsilon()) / std::mt19937::max();

	// compile time check!
	static_assert(std::mt19937::min() == 0, "hg_genrand_sp: std::mt19937::min() != 0");
	static_assert((float)std::mt19937::max() * fac < 1.0f, "hg_genrand_sp: < 1 failed");
	// Some codes have the following constant hard wired.
	static_assert((float)UINT32_MAX * 2.32830616e-10f < 1.0f, 
		      "hg_genrand_sp: check constant");	
	static_assert((float)UINT32_MAX * 2.32830616e-10f > 1.0f - 2.0f*std::numeric_limits<float>::epsilon(), 
		      "hg_genrand_sp: check constant");	

	auto y = (*rng)();
	*rn = (float)y * fac; /* reals: [0,1)-interval */
    }
}
