#include <bl_random_c.H>
#include <fstream>
#include <iostream>
#include <vector>

namespace 
{
    std::uint_fast32_t bl_rng_parallel_seed (int s, int rank, int nprocs)
    {
	std::seed_seq seq{s, 19937};
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

