#ifndef UTILS_H
#define UTILS_H

#include "particle.h"
#include <cmath>
#include <ctime>
#include <string>

double distance(Particle &particle1, Particle &particle2);

double distance2(Particle &particle1, Particle &particle2);

double calc_velocity(Particle& particle);

std::string get_current_time();

#endif