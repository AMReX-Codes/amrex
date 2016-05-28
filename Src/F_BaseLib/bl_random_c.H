#ifndef _BL_RANDOM_C_H_
#define _BL_RANDOM_C_H_

#include <random>

class BLRngUniformReal
{
public:
    BLRngUniformReal () {}
    BLRngUniformReal (std::uint_fast32_t s, double a, double b);
    ~BLRngUniformReal () {}
    double operator() ();
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::mt19937 m_eng;
    std::uniform_real_distribution<double> m_dist;
};

class BLRngNormal
{
public:
    BLRngNormal () {}
    BLRngNormal (std::uint_fast32_t s, double mean, double stddev);
    ~BLRngNormal () {}
    double operator() ();
    void save (const char* name) const;
    void restore (const char* name);
private:
    std::mt19937 m_eng;
    std::normal_distribution<double> m_dist;
};

class BLRngPoisson
{
public:
    BLRngPoisson () {}
    BLRngPoisson (std::uint_fast32_t s, double mean);
    ~BLRngPoisson () {}
    int operator() ();
    void save (const char* name) const;
    void restore (const char* name);
    void change_distribution (double mean);
private:
    std::mt19937 m_eng;
    std::poisson_distribution<int> m_dist;
};

class BLRngBinomial
{
public:
    BLRngBinomial () {}
    BLRngBinomial (std::uint_fast32_t s, int t, double p);
    ~BLRngBinomial () {}
    int operator() ();
    void save (const char* name) const;
    void restore (const char* name);
    void change_distribution (int t, double p);
private:
    std::mt19937 m_eng;
    std::binomial_distribution<int> m_dist;
};

#endif
