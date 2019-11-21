#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include "particle.hpp"
#include "neighborlist.hpp"
#include "unit_conversion.hpp"
#include "radial_distribution_function.hpp"
#include "property.hpp"
#include "utils.hpp"

using namespace std;

class Ensemble {

friend class Neighborlist;
friend class Rdf;
friend class Property;

public:
    // box(angstrom), temp(K), sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
    Ensemble(const unsigned _particle_number, double sigma, double epsilon, double mass, double init_temp, double set_temp, double time_interval, 
    double equilibration_time, double total_time, double box, string _output_path);

    // close file stream
    ~Ensemble();

    // return object of particle of that index
    inline Particle& operator[] (const int index);

    // calculation acceleration between a pair of particles
    inline void calc_acceleration(Particle& particle1, Particle& particle2);

    // correct the instaneous temperature by rescaling
    inline void rescale_temperature();

    inline void Andersen_thermostat(double collision_frequency);

    // main iteration
    inline void iteration();

    inline void energy_output(unsigned long i, ofstream& fout);
    inline void particle_movement_output(unsigned long i, Particle& particle, ofstream& fout);
    inline void temperature_output(unsigned long i, ofstream& fout);
    inline void msd_output(unsigned long i, double _MSD, ofstream& fout);

private:
    Unit_conversion unit;                           // Initialize unit conversion by sigma, epsilon, mass
    // all variable below are in reduced unit
    const double INIT_TEMP;                         // initial temperature before equilibration
    const double SET_TEMP;                          // set temperature after equilibration
    double TEMP;                                    // instaneous temperature
    const double BOX;                               // box dimension size
    const double TIME_INTERVAL;                     // machine time interval
    const double EQUILIBRATION_TIME;                // euquilibration time
    const double TOTAL_TIME;                        // total time
    const unsigned long EQUILIBRATION_ITERATION;    // iteration cycles of equilibration
    const unsigned long ITERATION;                  // total iteration cycles
    const unsigned long SAMPLE_RATE;                // sample rate defined as (iteration / 1000) such that the result contains 1000 points
    float ITERATION_PERCENTAGE;                      // percentage of main iteration
    unsigned particle_number;                       // particle number
    vector<Particle> ensemble;                      // main container of the particle ensemble
    const double rcut;                              // cutoff distance defined as 2.5 (reduced unit)
    const double ecut;                              // cutoff potential energy, calculated from 4.0 * (1 / pow(rcut, 12) - 1/ pow(rcut, 6))
    const double rlist2;                            // square of distance threshold of neighborlist defined as 3.5^2 (reduced unit)
    Neighborlist nlist;                             // object of neighborlist
    bool need_update_nlist;                         // whether or not to update the neighborlist
    Rdf rdf;                                        // object of radial distribution function
    Property property;                              // object of mean squared displacement calculation
    double ensemble_potential;                      // potential energy of ensemble
    double ensemble_kinetic;                        // kinetic energy of ensemble
    ofstream ensemble_out;                          // output file stream of energy
    ofstream particle_out;                          // output file stream of trajectory of selected particle
    ofstream temperature_out;                       // output file stream of temperature
    ofstream msd_out;                               // output file stream of mean square displacement
};


// box(angstrom), temp(K), sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
Ensemble::Ensemble(const unsigned _particle_number, double sigma, double epsilon, double mass, double init_temp, 
    double set_temp, double time_interval, double equilibration_time, double total_time, double box, string _output_path): 

    particle_number(_particle_number), unit(sigma, epsilon, mass), 
    INIT_TEMP(unit.reduced_temperature(init_temp)),
    SET_TEMP(unit.reduced_temperature(set_temp)), 
    TEMP(0),
    BOX(unit.reduced_distance(box * 1e-10)), 
    TIME_INTERVAL(unit.reduced_time(time_interval * 1e-15)), 
    EQUILIBRATION_TIME(unit.reduced_time(equilibration_time * 1e-9)),
    TOTAL_TIME(unit.reduced_time(total_time * 1e-9)), 
    EQUILIBRATION_ITERATION(EQUILIBRATION_TIME / TIME_INTERVAL),
    ITERATION(TOTAL_TIME / TIME_INTERVAL),
    SAMPLE_RATE(ITERATION / 1000),
    ensemble(_particle_number, Particle(TIME_INTERVAL)),
    rcut(2.5), ecut(-0.016316891), rlist2(12.25), 
    nlist(ensemble, BOX, rlist2), 
    rdf(1000, BOX, _output_path), 
    property((ITERATION - EQUILIBRATION_ITERATION) / SAMPLE_RATE + 1, unit.real_time(TIME_INTERVAL), _output_path),
    ensemble_out(_output_path + "/energy.csv"),
    particle_out(_output_path + "/particle.csv"),
    temperature_out(_output_path + "/temperature.csv"),
    msd_out(_output_path + "/msd.csv")
    {
        cout << "[MD LOG] " << get_current_time() << "\tMachine time interval: " << TIME_INTERVAL << endl;
        cout << "[MD LOG] " << get_current_time() << "\tEquilibration iteration: " << EQUILIBRATION_ITERATION << endl;
        cout << "[MD LOG] " << get_current_time() << "\tIteration: " << ITERATION << endl;
        cout << "[MD LOG] " << get_current_time() << "\tSample rate: " << SAMPLE_RATE << endl;
        cout << "[MD LOG] " << get_current_time() << "\tMachine box:" << BOX << endl;
        cout << "[MD LOG] " << get_current_time() << "\tEnsemble energy data output to \"" + _output_path + "/energy.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tParticle trajectory data output to \"" + _output_path + "/particle.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tTemperature data output to \"" + _output_path + "/temperature.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tDiffusion data output to \"" + _output_path + "/msd.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tDiffusion data output to \"" + _output_path + "/velocity_autocorr.csv\" ..." << endl;

        // initialize the postion and velocity of particles
        default_random_engine random_generator;
        uniform_real_distribution<double> displacement(0.0, 1.0);  //distribution generator
        normal_distribution<double> norm_dis(0.0, 1.0);
        double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0);
        double sumv2(0.0);
        int i = 0; // for lattice pos
        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            // particle->lattice_pos(i + 1);
            ++i;
            particle->pos_x += 0.01 * (displacement(random_generator) - 0.5); // cout << "Initial position: " << particle->pos_x << "\t";
            particle->pos_y += 0.01 * (displacement(random_generator) - 0.5); // cout << particle->pos_y << "\t";
            particle->pos_z += 0.01 * (displacement(random_generator) - 0.5); // cout << particle->pos_z << endl;
            particle->v_x = norm_dis(random_generator);
            particle->v_y = norm_dis(random_generator);
            particle->v_z = norm_dis(random_generator);

            sumv_x += particle->v_x;
            sumv_y += particle->v_y;
            sumv_z += particle->v_z;
            sumv2 += sumv_x * sumv_x + sumv_y * sumv_y + sumv_z * sumv_z;
        }
        sumv_x /= particle_number;
        sumv_y /= particle_number;
        sumv_z /= particle_number;
        sumv2 /= particle_number;
        double fs = sqrt(3 * INIT_TEMP / sumv2);
        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->v_x = (particle->v_x - sumv_x) * fs; // cout << "Initial velocity: " << particle->v_x << "\t";
            particle->v_y = (particle->v_y - sumv_y) * fs; // cout << particle->v_y << "\t";
            particle->v_z = (particle->v_z - sumv_z) * fs; // cout << particle->v_z << endl;
        }

        // Initialize neighborlist
        nlist.update_neighbor_list(ensemble);
        need_update_nlist = false;

        // Initialize acceleartion in step A
        for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
            for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                calc_acceleration(*particle1, *particle2);
                particle1->a_x_A = particle1->a_x_B;
                particle1->a_y_A = particle1->a_y_B;
                particle1->a_z_A = particle1->a_z_B;
                particle2->a_x_A = particle2->a_x_B;
                particle2->a_y_A = particle2->a_y_B;
                particle2->a_z_A = particle2->a_z_B;
            }
        }

        // execute movement and initialize a_B
        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->movement();
            particle->velocity();
            particle->a_x_A = particle->a_x_B;
            particle->a_y_A = particle->a_y_B;
            particle->a_z_A = particle->a_z_B;
            particle->a_x_B = 0;
            particle->a_y_B = 0;
            particle->a_z_B = 0;
            particle->potential_value = 0;
        }
}


Ensemble::~Ensemble() {
    ensemble_out.close();
    particle_out.close();
    temperature_out.close();
    msd_out.close();
    cout << "[MD LOG] " << get_current_time() << "\tOutput file saved" << endl;
}


// use ensemble[index] to call particle's method
Particle& Ensemble::operator[] (const int index) {
    return ensemble[index];
}


void Ensemble::calc_acceleration(Particle& particle1, Particle& particle2) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 <= rcut * rcut) {
        double r2i = 1 / r2;
        double r6i = pow(r2i, 3);
        double force_x = 48.0 * r6i * r2i * (r6i - 0.5) * dx;
        double force_y = 48.0 * r6i * r2i * (r6i - 0.5) * dy;
        double force_z = 48.0 * r6i * r2i * (r6i - 0.5) * dz;
        particle1.a_x_B += force_x;
        particle1.a_y_B += force_y;
        particle1.a_z_B += force_z;
        particle2.a_x_B -= force_x;
        particle2.a_y_B -= force_y;
        particle2.a_z_B -= force_z;
        particle1.potential_value += (4.0 * r6i * ( r6i - 1 ) - ecut);
    }
    if (r2 > nlist.rlist2) {
        need_update_nlist = true;
    }
}


void Ensemble::rescale_temperature() {
    double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0);
    double sumv2(0.0);
    for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
        sumv_x += particle->v_x;
        sumv_y += particle->v_y;
        sumv_z += particle->v_z;
        sumv2 += sumv_x * sumv_x + sumv_y * sumv_y + sumv_z * sumv_z;
    }
    sumv_x /= particle_number;
    sumv_y /= particle_number;
    sumv_z /= particle_number;
    sumv2 /= particle_number;
    TEMP = sumv2 / 3;
    double fs = sqrt(SET_TEMP / TEMP);
    for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
        particle->v_x = (particle->v_x - sumv_x) * fs;
        particle->v_y = (particle->v_y - sumv_y) * fs;
        particle->v_z = (particle->v_z - sumv_z) * fs;
    }
}


void Ensemble::Andersen_thermostat(double collision_frequency) {
    double sigma = sqrt(SET_TEMP);
    double scale_factor;
    default_random_engine random_generator;
    normal_distribution<double> gauss(0, sigma);
    uniform_real_distribution<double> ranf(0.0, 1.0);
    for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
        if (ranf(random_generator) < collision_frequency) {
            scale_factor = gauss(random_generator) / calc_velocity(*it);
            it->v_x *= scale_factor;
            it->v_y *= scale_factor;
            it->v_z *= scale_factor;
            scale_factor = 0;
        }
    }
}


// void Ensemble::energy_output(unsigned long i, ofstream& fout) {
//     fout << i * unit.real_time(TIME_INTERVAL) << "    " << unit.real_energy(ensemble_potential) 
//         << "    " << unit.real_energy(ensemble_kinetic) << "    " 
//         << unit.real_energy(ensemble_potential + ensemble_kinetic) << endl;
// } 
void Ensemble::energy_output(unsigned long i, ofstream& fout) {
    fout << i << "    " << ensemble_potential << "    " << ensemble_kinetic << "    " 
        << ensemble_potential + ensemble_kinetic << endl;
}


// void Ensemble::particle_movement_output(unsigned long i, Particle& particle, ofstream& fout) {
//     fout << i * unit.real_time(TIME_INTERVAL) << "    " << unit.real_distance(particle.pos_x) << "    " 
//         << unit.real_distance(particle.pos_y) << "    " << unit.real_distance(particle.pos_z) << "    " 
//         << unit.real_velocity(particle.v_x) << "    " << unit.real_velocity(particle.v_y) << "    " 
//         << unit.real_velocity(particle.v_z) << "    " << unit.real_acceleration(particle.a_x_A) << "    " 
//         << unit.real_acceleration(particle.a_y_A) << "    " << unit.real_acceleration(particle.a_z_A) << "    " 
//         << unit.real_energy(particle.potential_value) << "    " << unit.real_energy(particle.kinetic_value) 
//         << "    " << unit.real_energy(particle.potential_value + particle.kinetic_value) << endl;
// }

void Ensemble::particle_movement_output(unsigned long i, Particle& particle, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << particle.pos_x << "    " << particle.pos_y << "    " << particle.pos_z << "    " 
        << particle.v_x << "    " << particle.v_y << "    " << particle.v_z << "    " << particle.a_x_A << "    " 
        << particle.a_y_A << "    " << particle.a_z_A << "    " << particle.potential_value << "    " << particle.kinetic_value 
        << "    " << particle.potential_value + particle.kinetic_value << endl;
}


void Ensemble::temperature_output(unsigned long i, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << TEMP << endl;
}

void Ensemble::msd_output(unsigned long i, double _MSD, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << pow(sqrt(_MSD), 2) << endl;
}

void Ensemble::iteration() {
    unsigned long i = 0;
    while (i <= ITERATION) {
        // Initialize ensemble energy
        ensemble_kinetic = 0;
        ensemble_potential = 0;

        // calculate acceleration of step B in neighbor list
        for (int i = 0; i < nlist.nlist.size(); ++i) {
            for (auto j = nlist.nlist[i].begin(); j != nlist.nlist[i].end(); ++j) {
                calc_acceleration(ensemble[i], ensemble[*j]);
            }
        }

        if (need_update_nlist == true) {
            nlist.update_neighbor_list(ensemble);
            need_update_nlist = false;
        }

        // initialize summary of velocity
        double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0);
        double sumv2(0.0);
        // calculate velocity of step B
        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            ensemble_potential += particle->potential_value;

            // get the instaneous temperature
            sumv_x += particle->v_x;
            sumv_y += particle->v_y;
            sumv_z += particle->v_z;
            sumv2 += sumv_x * sumv_x + sumv_y * sumv_y + sumv_z * sumv_z;

            // execute x and v propagation
            particle->movement();
            particle->velocity();

            // record a_A and initialize a_B
            particle->a_x_A = particle->a_x_B;
            particle->a_y_A = particle->a_y_B;
            particle->a_z_A = particle->a_z_B;
            particle->a_x_B = 0;
            particle->a_y_B = 0;
            particle->a_z_B = 0;
            particle->potential_value = 0;
        }
        ensemble_kinetic = 0.5 * sumv2;
        TEMP = sumv2 / (3 * particle_number);

        // after equilibration iteration of initial equilibrition, control temperature by rescaling
        if (i == EQUILIBRATION_ITERATION) {
            rescale_temperature();
            property.initalize(ensemble);
        }

        if (i >= EQUILIBRATION_ITERATION) {
            // Andersen_thermostat(1e-4);
            if (i % SAMPLE_RATE == 0) {
                particle_movement_output(i, ensemble[1], particle_out);
                energy_output(i, ensemble_out);
                rdf.sample(ensemble);
            }
            unsigned long VACF_sample_iteration = 5000;
            if (EQUILIBRATION_ITERATION + 10000 + VACF_sample_iteration >= ITERATION) {
                cerr << "[MD ERR]\t" << get_current_time() << "\tIncorrect VACF sample iteration" << endl;
                break;
            }
            if (i > EQUILIBRATION_ITERATION + 10000 && i <= EQUILIBRATION_ITERATION + 10000 + VACF_sample_iteration) {
                msd_output(i, property.calc_mean_square_particle_displacement(ensemble), msd_out);
                property.sample_velocity_autocorrelation(ensemble);
            }
        }

        if (i % SAMPLE_RATE == 0) {
            ITERATION_PERCENTAGE = ((float) i / (float) ITERATION) * 100;
            cout << "\r[MD LOG] " << get_current_time() << "\t" << ITERATION_PERCENTAGE << "\% completed " << flush;
            temperature_output(i, temperature_out);
        }
        ++i;
    }
    cout << endl;   // output a new line for the progess log
    rdf.normalize(particle_number);
    rdf.output();
    property.calc_velocity_autocorrelation();
    property.velocity_autocorr_output();
}


#endif
