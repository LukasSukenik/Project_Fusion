#ifndef RNG_H
#define RNG_H

#include <random>

typedef long double myFloat;

std::mt19937_64 rng;
std::uniform_real_distribution<double> unif;

#endif // RNG_H
