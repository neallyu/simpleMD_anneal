#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

#include <cmath>

class Unit_conversion {

public:
    // receive format: sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
    Unit_conversion(double _sigma, double _epsilon, double _mass): 
        SIGMA(_sigma * 1e-10), EPSILON(_epsilon * 1000 / NA), MASS(_mass / (1000 * NA)) { }

    double real_distance(double _distance) {
        return _distance * SIGMA;
    }

    double reduced_distance(double _distance) {
        return _distance / SIGMA;
    }

    double real_velocity(double _velocity) {
        return sqrt(EPSILON / MASS) * _velocity;
    }

    double real_acceleration(double _acceleration) {
        return _acceleration * EPSILON / (MASS * SIGMA);
    }

    double real_time(double _time) {
        return SIGMA * sqrt(MASS / EPSILON) * _time;
    }

    double reduced_time(double _time) {
        return _time * sqrt(EPSILON / MASS) / SIGMA;
    }

    double real_temperature(double _temperature) {
        return EPSILON * _temperature / kb;
    }

    double reduced_temperature(double _temperature) {
        return kb *_temperature / EPSILON;
    }

    double real_pressure(double _pressure) {
        return EPSILON * _pressure / pow(SIGMA, 3);
    }

    double real_density(double _density) {
        return _density / pow(SIGMA, 3);
    }

    double real_energy(double _energy) {
        return _energy * EPSILON;
    }

private:
    const double kb = 1.380649e-23; // bolzmann constant (J/K)
    const double NA = 6.02214076e23; // Avogadro constant

    const double SIGMA;
    const double EPSILON;
    const double MASS;
};

#endif