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
    ofs << rng;
}

template <typename T>
void
BLRng_restore (T& rng, const char* name)
{
    std::ifstream ifs(name);
    if (ifs.good()) {
	ifs >> rng;
    } else {
	std::cerr << "bl_rng: faied to open " << name << std::endl;
	backtrace_handler(6);
    }
}

template <typename Return_t, typename Distro_t, typename Engine_t>
Return_t 
BLRng_get (Distro_t& distro, Engine_t& eng)
{
    return distro(eng);
}

#endif
