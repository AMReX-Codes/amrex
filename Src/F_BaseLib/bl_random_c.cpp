#include <bl_random_c.H>
#include <fstream>
#include <iostream>
#include <vector>
#include <set>
#include <limits>

extern "C" { void backtrace_handler(int); }

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
}

// 
// engine
//
extern "C"
{
    void bl_rng_new_engine_c (BLRngEngine*& eng, int s, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	eng = new BLRngEngine(seed);
    }

    void bl_rng_delete_engine_c (BLRngEngine* eng)
    {
	delete eng;
    }

    void bl_rng_save_engine_c (const BLRngEngine* eng, const char* name)
    {
	eng->save(name);
    }

    void bl_rng_restore_engine_c (BLRngEngine*& eng, const char* name)
    {
	eng = new BLRngEngine();
	eng->restore(name);
    }
}

BLRngEngine::BLRngEngine (std::uint_fast32_t s)
    : m_eng(s)
{
    m_eng.discard(1000000);  // warm up
}

void
BLRngEngine::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_eng;
}

void
BLRngEngine::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_eng;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

//
// uniform real distribution
//
extern "C"
{
    void bl_rng_new_uniform_real_c (BLRngUniformReal*& rng, double a, double b)
    {
	rng = new BLRngUniformReal(a,b);
    }
    //
    void bl_rng_delete_uniform_real_c (BLRngUniformReal* rng)
    {
	delete rng;
    }
    //
    double bl_rng_get_uniform_real_c (BLRngUniformReal* rng, BLRngEngine* eng)
    {
	return (*rng)(*eng);
    }
    //
    void bl_rng_save_uniform_real_c (const BLRngUniformReal* rng, const char* name)
    {
	rng->save(name);
    }
    //
    void bl_rng_restore_uniform_real_c (BLRngUniformReal*& rng, const char* name)
    {
	rng = new BLRngUniformReal();
	rng->restore(name);
    }
}

BLRngUniformReal::BLRngUniformReal (double a, double b)
    : m_dist(a,b)
{ }

double
BLRngUniformReal::operator() (BLRngEngine& eng)
{
    return m_dist(eng.get());
}

void
BLRngUniformReal::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_dist;
}

void
BLRngUniformReal::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_dist;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

//
// normal distribution
//
extern "C"
{
    void bl_rng_new_normal_c (BLRngNormal*& rng, double mean, double stddev)
    {
	rng = new BLRngNormal(mean,stddev);
    }
    //
    void bl_rng_delete_normal_c (BLRngNormal* rng)
    {
	delete rng;
    }
    //
    double bl_rng_get_normal_c (BLRngNormal* rng, BLRngEngine* eng)
    {
	return (*rng)(*eng);
    }
    //
    void bl_rng_save_normal_c (const BLRngNormal* rng, const char* name)
    {
	rng->save(name);
    }
    //
    void bl_rng_restore_normal_c (BLRngNormal*& rng, const char* name)
    {
	rng = new BLRngNormal();
	rng->restore(name);
    }
}

BLRngNormal::BLRngNormal (double mean, double stddev)
    : m_dist(mean,stddev)
{ }

double
BLRngNormal::operator() (BLRngEngine& eng)
{
    return m_dist(eng.get());
}

void
BLRngNormal::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_dist;
}

void
BLRngNormal::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_dist;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

//
// poisson distribution
//
extern "C"
{
    void bl_rng_new_poisson_c (BLRngPoisson*& rng, double mean)
    {
	rng = new BLRngPoisson(mean);
    }
    //
    void bl_rng_delete_poisson_c (BLRngPoisson* rng)
    {
	delete rng;
    }
    //
    int bl_rng_get_poisson_c (BLRngPoisson* rng, BLRngEngine* eng)
    {
	return (*rng)(*eng);
    }
    //
    void bl_rng_save_poisson_c (const BLRngPoisson* rng, const char* name)
    {
	rng->save(name);
    }
    //
    void bl_rng_restore_poisson_c (BLRngPoisson*& rng, const char* name)
    {
	rng = new BLRngPoisson();
	rng->restore(name);
    }

    void bl_rng_change_poisson_c (BLRngPoisson* rng, double mean)
    {
	rng->change_distribution(mean);
    }
}

BLRngPoisson::BLRngPoisson (double mean)
    : m_dist(mean)
{ }

int
BLRngPoisson::operator() (BLRngEngine& eng)
{
    return m_dist(eng.get());
}

void
BLRngPoisson::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_dist;
}

void
BLRngPoisson::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_dist;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

void
BLRngPoisson::change_distribution (double mean)
{
    if (mean != m_dist.mean())
	m_dist = std::poisson_distribution<int>(mean);
}

//
// binomial distribution
//
extern "C"
{
    void bl_rng_new_binomial_c (BLRngBinomial*& rng, int t, double p)
    {
	rng = new BLRngBinomial(t,p);
    }
    //
    void bl_rng_delete_binomial_c (BLRngBinomial* rng)
    {
	delete rng;
    }
    //
    int bl_rng_get_binomial_c (BLRngBinomial* rng, BLRngEngine* eng)
    {
	return (*rng)(*eng);
    }
    //
    void bl_rng_save_binomial_c (const BLRngBinomial* rng, const char* name)
    {
	rng->save(name);
    }
    //
    void bl_rng_restore_binomial_c (BLRngBinomial*& rng, const char* name)
    {
	rng = new BLRngBinomial();
	rng->restore(name);
    }
    //
    void bl_rng_change_binomial_c (BLRngBinomial* rng,
				   int t, double p)
    {
	rng->change_distribution(t,p);
    }
}

BLRngBinomial::BLRngBinomial (int t, double p)
    : m_dist(t,p)
{ }

int
BLRngBinomial::operator() (BLRngEngine& eng)
{
    return m_dist(eng.get());
}

void
BLRngBinomial::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_dist;
}

void
BLRngBinomial::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_dist;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

void
BLRngBinomial::change_distribution (int t, double p)
{
    if (t != m_dist.t() || p != m_dist.p())
	m_dist = std::binomial_distribution<int>(t,p);
}
