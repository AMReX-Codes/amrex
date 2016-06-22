#ifndef _BL_RANDOM_C_H_
#define _BL_RANDOM_C_H_

#include <random>

class BLRngEngine
{
public:
    BLRngEngine () = default;
    BLRngEngine (std::uint_fast32_t s);
    ~BLRngEngine () {}
    std::mt19937& get () { return m_eng; }
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::mt19937 m_eng;
};

class BLRngUniformReal
{
public:
    BLRngUniformReal () = default;
    BLRngUniformReal (double a, double b);
    ~BLRngUniformReal () {}
    double operator() (BLRngEngine& eng);
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::uniform_real_distribution<double> m_dist;
};

class BLRngNormal
{
public:
    BLRngNormal () = default;
    BLRngNormal (double mean, double stddev);
    ~BLRngNormal () {}
    double operator() (BLRngEngine& eng);
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::normal_distribution<double> m_dist;
};

class BLRngPoisson
{
public:
    BLRngPoisson () = default;
    BLRngPoisson (double mean);
    ~BLRngPoisson () {}
    int operator() (BLRngEngine& eng);
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::poisson_distribution<int> m_dist;
};

class BLRngBinomial
{
public:
    BLRngBinomial () = default;
    BLRngBinomial (int t, double p);
    ~BLRngBinomial () {}
    int operator() (BLRngEngine& eng);
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::binomial_distribution<int> m_dist;
};

#endif
