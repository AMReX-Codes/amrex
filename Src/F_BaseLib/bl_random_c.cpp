#include <bl_random_c.H>
#include <fstream>
#include <iostream>
#include <vector>

namespace 
{
    std::uint_fast32_t bl_rng_parallel_seed (int s, int rank, int nprocs)
    {
	std::seed_seq seq{s, s+4208, s-76, s-17, s+19937, s+1};
	std::vector<std::uint32_t> seeds(nprocs);
	seq.generate(seeds.begin(), seeds.end());
	return seeds[rank];
    }
}

extern "C" { void backtrace_handler(int); }

//
// uniform real distribution
//
extern "C"
{
    void bl_rng_new_uniform_real_c (BLRngUniformReal*& rng,
				    int s, double a, double b, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	rng = new BLRngUniformReal(seed,a,b);
    }
    //
    void bl_rng_delete_uniform_real_c (BLRngUniformReal* rng)
    {
	delete rng;
    }
    //
    double bl_rng_get_uniform_real_c (BLRngUniformReal* rng)
    {
	return (*rng)();
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

BLRngUniformReal::BLRngUniformReal (std::uint_fast32_t s, double a, double b)
    : m_eng(s), m_dist(a,b)
{
    // improve quality of poor seeds
    m_eng.discard(1000000);
}

double
BLRngUniformReal::operator() ()
{
    return m_dist(m_eng);
}

void
BLRngUniformReal::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_eng << "\n" << m_dist << "\n";
}

void
BLRngUniformReal::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_eng >> m_dist;
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
    void bl_rng_new_normal_c (BLRngNormal*& rng,
			      int s, double mean, double stddev, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	rng = new BLRngNormal(seed,mean,stddev);
    }
    //
    void bl_rng_delete_normal_c (BLRngNormal* rng)
    {
	delete rng;
    }
    //
    double bl_rng_get_normal_c (BLRngNormal* rng)
    {
	return (*rng)();
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

BLRngNormal::BLRngNormal (std::uint_fast32_t s, double mean, double stddev)
    : m_eng(s), m_dist(mean,stddev)
{
    // improve quality of poor seeds
    m_eng.discard(1000000);
}

double
BLRngNormal::operator() ()
{
    return m_dist(m_eng);
}

void
BLRngNormal::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_eng << "\n" << m_dist << "\n";
}

void
BLRngNormal::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_eng >> m_dist;
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
    void bl_rng_new_poisson_c (BLRngPoisson*& rng,
			      int s, double mean, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	rng = new BLRngPoisson(seed,mean);
    }
    //
    void bl_rng_delete_poisson_c (BLRngPoisson* rng)
    {
	delete rng;
    }
    //
    int bl_rng_get_poisson_c (BLRngPoisson* rng)
    {
	return (*rng)();
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

    void bl_rng_change_poisson_c (BLRngPoisson* rng,
				  double mean)
    {
	rng->change_distribution(mean);
    }
}

BLRngPoisson::BLRngPoisson (std::uint_fast32_t s, double mean)
    : m_eng(s), m_dist(mean)
{
    // improve quality of poor seeds
    m_eng.discard(1000000);
}

int
BLRngPoisson::operator() ()
{
    return m_dist(m_eng);
}

void
BLRngPoisson::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_eng << "\n" << m_dist << "\n";
}

void
BLRngPoisson::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_eng >> m_dist;
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
    void bl_rng_new_binomial_c (BLRngBinomial*& rng,
				int s, int t, double p, int rank, int nprocs)
    {
	std::uint_fast32_t seed = bl_rng_parallel_seed(s, rank, nprocs);
	rng = new BLRngBinomial(seed,t,p);
    }
    //
    void bl_rng_delete_binomial_c (BLRngBinomial* rng)
    {
	delete rng;
    }
    //
    int bl_rng_get_binomial_c (BLRngBinomial* rng)
    {
	return (*rng)();
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

BLRngBinomial::BLRngBinomial (std::uint_fast32_t s, int t, double p)
    : m_eng(s), m_dist(t,p)
{
    // improve quality of poor seeds
    m_eng.discard(1000000);
}

int
BLRngBinomial::operator() ()
{
    return m_dist(m_eng);
}

void
BLRngBinomial::save (const char* name) const
{
    std::ofstream ofs(name);
    ofs << m_eng << "\n" << m_dist << "\n";
}

void
BLRngBinomial::restore (const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> m_eng >> m_dist;
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
