#ifndef UTILS_H
#define UTILS_H

#include "particle.hpp"
#include <cmath>
#include <ctime>
#include <string>

double distance(Particle &particle1, Particle &particle2, double BOX) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    return sqrt(dx * dx + dy * dy + dz * dz);
}


double distance2(Particle &particle1, Particle &particle2, double BOX) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    return dx * dx + dy * dy + dz * dz;
}


double calc_velocity(Particle& particle) {
    return sqrt(pow(particle.v_x, 2) + pow(particle.v_y, 2) + pow(particle.v_z, 2));
}


std::string get_current_time() {
    time_t now = std::time(0);
    tm *ltm = std::localtime(&now);
    std::string current_time;
    char space[] = " ";
    char bar[] = "-";
    char colon[] = ":";
    char zero[] = "0";

    current_time.append(std::to_string(1900 + ltm->tm_year));
    current_time.push_back(*bar);
    current_time.append(std::to_string(1 + ltm->tm_mon));
    current_time.push_back(*bar);
    current_time.append(std::to_string(ltm->tm_mday));
    current_time.push_back(*space);
    current_time.append(std::to_string(ltm->tm_hour));
    current_time.push_back(*colon);
    if (ltm->tm_min < 10) {
        current_time.push_back(*zero);
        current_time.append(std::to_string(ltm->tm_min));
    } else {
        current_time.append(std::to_string(ltm->tm_min));
    }
    current_time.push_back(*colon);
    if (ltm->tm_sec < 10) {
        current_time.push_back(*zero);
        current_time.append(std::to_string(ltm->tm_sec));
    } else {
        current_time.append(std::to_string(ltm->tm_sec));
    }

    return current_time;
}

#endif