#ifndef _BL_RANDOM_C_H_
#define _BL_RANDOM_C_H_

#include <random>
#include <fstream>
#include <iostream>

extern "C" { void backtrace_handler(int); }

template <typename T>
void
BLRng_save (const T& rng, const char* name)
{
    std::ofstream ofs(name);
    ofs << rng.m_data;
}

template <typename T>
void
BLRng_restore (T& rng, const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> rng.m_data;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

class BLRngEngine
{
public:
    BLRngEngine () = default;
    BLRngEngine (std::uint_fast32_t s);
    ~BLRngEngine () {}
    std::mt19937& get () { return m_data; }
    template <typename T> friend void BLRng_save(const T& rng, const char* name);
    template <typename T> friend void BLRng_restore(T& rng, const char* name);
private:
    std::mt19937 m_data;
};

class BLRngUniformReal
{
public:
    BLRngUniformReal () = default;
    BLRngUniformReal (double a, double b);
    ~BLRngUniformReal () {}
    double operator() (BLRngEngine& eng);
    template <typename T> friend void BLRng_save(const T& rng, const char* name);
    template <typename T> friend void BLRng_restore(T& rng, const char* name);
private:
    std::uniform_real_distribution<double> m_data;
};

class BLRngNormal
{
public:
    BLRngNormal () = default;
    BLRngNormal (double mean, double stddev);
    ~BLRngNormal () {}
    double operator() (BLRngEngine& eng);
    template <typename T> friend void BLRng_save(const T& rng, const char* name);
    template <typename T> friend void BLRng_restore(T& rng, const char* name);
private:
    std::normal_distribution<double> m_data;
};

class BLRngPoisson
{
public:
    BLRngPoisson () = default;
    BLRngPoisson (double mean);
    ~BLRngPoisson () {}
    int operator() (BLRngEngine& eng);
    template <typename T> friend void BLRng_save(const T& rng, const char* name);
    template <typename T> friend void BLRng_restore(T& rng, const char* name);
private:
    std::poisson_distribution<int> m_data;
};

class BLRngBinomial
{
public:
    BLRngBinomial () = default;
    BLRngBinomial (int t, double p);
    ~BLRngBinomial () {}
    int operator() (BLRngEngine& eng);
    template <typename T> friend void BLRng_save(const T& rng, const char* name);
    template <typename T> friend void BLRng_restore(T& rng, const char* name);
private:
    std::binomial_distribution<int> m_data;
};

#endif
