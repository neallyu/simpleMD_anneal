#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <omp.h>
#include "particle.h"
#include "utils.h"

using namespace std;

class Ensemble {

public:
    // reduced unit
    Ensemble(const unsigned _particle_number, double _box, double init_temp, double set_temp, double time_interval, 
    double equilibration_time, double total_time, string _output_path);

    // close file stream
    ~Ensemble();

    // return object of particle of that index
    inline Particle& operator[] (const int index);

    // lattice position
    inline void lattice_pos_box();
    inline void lattice_pos_sphere();

    // calculation acceleration between a pair of particles
    inline void calc_acceleration(Particle& particle1, Particle& particle2);

    inline void box_boundary(Particle& particle);

    inline void sphere_boundary(Particle& particle);

    inline double temperature_decreasing_curve(unsigned long i);

    // correct the instaneous temperature by rescaling
    inline void rescale_temperature(double target_TEMP);

    inline void Andersen_thermostat(double targetTEMP, double collision_frequency);

    inline void recenter();
    inline void recenter_sphere();

    // main iteration
    void iteration();

    inline void energy_output(unsigned long i, ofstream& fout);
    inline void particle_movement_output(unsigned long i, Particle& particle, ofstream& fout);
    inline void temperature_output(unsigned long i, ofstream& fout);
    inline void coordinates_output(ofstream& fout_cif, ofstream& fout2_rawdata);

private:
    // all variable below are in reduced unit
    const double INIT_TEMP;                         // initial temperature before equilibration
    const double SET_TEMP;                          // set temperature after equilibration
    double TEMP;                                    // instaneous temperature
    unsigned particle_number;                       // particle number
    const double BOX;                               // box dimension size
    const double TIME_INTERVAL;                     // machine time interval
    const double EQUILIBRATION_TIME;                // euquilibration time
    const double TOTAL_TIME;                        // total time
    const unsigned long EQUILIBRATION_ITERATION;    // iteration cycles of equilibration
    const double EQUILIBRATION_ITERATION_2;         // square of EQUILIBRATION_ITERATION
    const unsigned long ITERATION;                  // total iteration cycles
    const unsigned long SAMPLE_RATE;                // sample rate defined as (iteration / 1000) such that the result contains 1000 points
    float ITERATION_PERCENTAGE;                      // percentage of main iteration
    vector<Particle> ensemble;                      // main container of the particle ensemble
    double ensemble_potential;                      // potential energy of ensemble
    double ensemble_kinetic;                        // kinetic energy of ensemble
    ofstream ensemble_out;                          // output file stream of energy
    ofstream particle_out;                          // output file stream of trajectory of selected particle
    ofstream temperature_out;                       // output file stream of temperature
    ofstream coordinates_cif_out;                   // output file stream of coordinates of particles
    ofstream coordinates_rawdata_out;               
};


#endif