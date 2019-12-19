#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include <omp.h>
#include "particle.hpp"
#include "utils.hpp"
#include "conjugated_gradient.hpp"

using namespace std;

class Ensemble {

friend class CG;

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
    inline void iteration();

    inline void energy_output(unsigned long i, ofstream& fout);
    inline void particle_movement_output(unsigned long i, Particle& particle, ofstream& fout);
    inline void temperature_output(unsigned long i, ofstream& fout);
    inline void coordinates_output(ofstream& fout);

    inline void inter_distance_output(unsigned long i, ofstream& fout);

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
    ofstream coordinates_out;                       // output file stream of coordinates of particles
    string outputPath;                              // store the output path
    ofstream inter_distance_out;                    // output file stream of inter particle distance summary
};

// reduced unit
Ensemble::Ensemble(const unsigned _particle_number, double _box, double init_temp, double set_temp, double time_interval, 
    double equilibration_time, double total_time, string _output_path): 
    particle_number(_particle_number),
    INIT_TEMP(init_temp),
    SET_TEMP(set_temp), 
    TEMP(0),
    BOX(_box), 
    TIME_INTERVAL(time_interval), 
    EQUILIBRATION_TIME(equilibration_time),
    TOTAL_TIME(total_time), 
    EQUILIBRATION_ITERATION(EQUILIBRATION_TIME / TIME_INTERVAL),
    EQUILIBRATION_ITERATION_2(pow(EQUILIBRATION_ITERATION, 2)),
    ITERATION(TOTAL_TIME / TIME_INTERVAL),
    SAMPLE_RATE(ITERATION / 1000),
    ensemble(_particle_number, Particle(TIME_INTERVAL)),
    ensemble_out(_output_path + "/energy.csv"),
    particle_out(_output_path + "/particle.csv"),
    temperature_out(_output_path + "/temperature.csv"),
    coordinates_out(_output_path + "/coordinates.cif"),
    inter_distance_out(_output_path + "/inter_distance.csv"),
    outputPath(_output_path)
    {
        cout << "[MD LOG] " << get_current_time() << "\tEquilibration iteration: " << EQUILIBRATION_ITERATION << endl;
        cout << "[MD LOG] " << get_current_time() << "\tIteration: " << ITERATION << endl;
        cout << "[MD LOG] " << get_current_time() << "\tSample rate: " << SAMPLE_RATE << endl;
        cout << "[MD LOG] " << get_current_time() << "\tBox: " << BOX << endl;
        cout << "[MD LOG] " << get_current_time() << "\tEnsemble energy data output to \"" + _output_path + "/energy.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tParticle trajectory data output to \"" + _output_path + "/particle.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tTemperature data output to \"" + _output_path + "/temperature.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tCoordinates data output to \"" + _output_path + "/coordinates.cif\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tCoordinates data output to \"" + _output_path + "/coordinates_after_CG.cif\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tCoordinates data output to \"" + _output_path + "/potential_gradient_norm.csv\" ..." << endl;

        // parallel
        omp_set_num_threads(32);

        lattice_pos_box();
        // lattice_pos_sphere();
        rescale_temperature(INIT_TEMP);

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
    coordinates_out.close();
    inter_distance_out.close();
    cout << "[MD LOG] " << get_current_time() << "\tOutput file saved" << endl;
}


// use ensemble[index] to call particle's method
Particle& Ensemble::operator[] (const int index) {
    return ensemble[index];
}


void Ensemble::lattice_pos_box() {

    default_random_engine random_generator;
    uniform_real_distribution<double> displacement(0, 1.0);  //distribution generator
    normal_distribution<double> norm_dis(0.0, 1.0);

    for (int n = 0; n < particle_number; ++n) {
        int i;
        for (i = 0; pow(i, 3) < n + 1; ++i);
        int j = n + 1 - (i - 1) * (i - 1) * (i - 1);
        if (j <= i * i) {
            ensemble[n].pos_x = i - 1;
            ensemble[n].pos_y = (j - 1) / i;
            ensemble[n].pos_z = (j - 1) % i;
        }
        else if (j <= i * (2 * i - 1)) {
            ensemble[n].pos_x = (j - i * i - 1) / i;
            ensemble[n].pos_y = i - 1;
            ensemble[n].pos_z = (j - i * i - 1) % i;
        }
        else {
            ensemble[n].pos_x = (j - i * (2 * i - 1) - 1) / (i - 1);
            ensemble[n].pos_y = (j - i * (2 * i - 1) - 1) % (i - 1);
            ensemble[n].pos_z = i - 1;
        }
        ensemble[n].pos_x += 0.01 * displacement(random_generator) ;
        ensemble[n].pos_y += 0.01 * displacement(random_generator) ;
        ensemble[n].pos_z += 0.01 * displacement(random_generator) ;

        // cout << "pos_x: " << ensemble[n].pos_x << "\tpos_y: " << ensemble[n].pos_y << "\tpos_z: " << ensemble[n].pos_z << endl;

        ensemble[n].v_x = norm_dis(random_generator);
        ensemble[n].v_y = norm_dis(random_generator);
        ensemble[n].v_z = norm_dis(random_generator);
    }
}


void Ensemble::lattice_pos_sphere() {
    lattice_pos_box();
    recenter_sphere();
}


void Ensemble::calc_acceleration(Particle& particle1, Particle& particle2) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    double r2 = dx * dx + dy * dy + dz * dz;
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
        particle1.potential_value += 4.0 * r6i * (r6i - 1);
        // particle2.potential_value += 4.0 * r6i * (r6i - 1);
}


void Ensemble::box_boundary(Particle &particle) {
    if (particle.pos_x < 0 && particle.v_x < 0) {
        particle.v_x *= -1;
    }
    if (particle.pos_y < 0 && particle.v_y < 0) {
        particle.v_y *= -1;
    }
    if (particle.pos_z < 0 && particle.v_z < 0) {
        particle.v_z *= -1;
    }
    if (particle.pos_x > BOX && particle.v_x > 0) {
        particle.v_x *= -1;
    }
    if (particle.pos_y > BOX && particle.v_y > 0) {
        particle.v_y *= -1;
    }
    if (particle.pos_z > BOX && particle.v_z > 0) {
        particle.v_z *= -1;
    }
}


void Ensemble::sphere_boundary(Particle& particle) {
    double r2 = pow(particle.pos_x, 2) + pow(particle.pos_y, 2) + pow(particle.pos_z, 2);
    double BOX_2 = BOX * BOX;
    if (r2 >= BOX_2) {
        vector<double> pos_vec({particle.pos_x, particle.pos_y, particle.pos_z});
        vector<double> velocity_vec({particle.v_x, particle.v_y, particle.v_z});

        //scale factor to make new y-axis of linear relavant vector to pos vec and velocity vec which is perpendicular to the new x-axis
        double a = -(pow(pos_vec[0], 2) + pow(pos_vec[1], 2) + pow(pos_vec[2], 2)) / 
            (pos_vec[0] * velocity_vec[0] + pos_vec[1] * velocity_vec[1] + pos_vec[2] * velocity_vec[2]);

        vector<double> y_vector(pos_vec);
        y_vector[0] *= a; y_vector[0] += velocity_vec[0];
        y_vector[1] *= a; y_vector[1] += velocity_vec[1];
        y_vector[2] *= a; y_vector[2] += velocity_vec[2];

        vector<double> z_vector(3, 0);
        z_vector[0] = pos_vec[1] * y_vector[2] - pos_vec[2] * y_vector[1];
        z_vector[1] = pos_vec[2] * y_vector[0] - pos_vec[0] * y_vector[2];
        z_vector[2] = pos_vec[0] * y_vector[1] - pos_vec[1] * y_vector[0];

        // linear transformation of velocity vector from original cooridinates to new cooridnates
        vector<double> new_velocity_vec(3, 0);
        new_velocity_vec[0] = pos_vec[0] * velocity_vec[0] + y_vector[0] * velocity_vec[1] + z_vector[0] * velocity_vec[2];
        new_velocity_vec[1] = pos_vec[1] * velocity_vec[0] + y_vector[1] * velocity_vec[1] + z_vector[1] * velocity_vec[2];
        new_velocity_vec[2] = pos_vec[2] * velocity_vec[0] + y_vector[2] * velocity_vec[1] + z_vector[2] * velocity_vec[2];

        // if the x-portion of the vector larger than zero, do mirror operation to the x-portion of transformed velocity vector
        if (new_velocity_vec[0] > 0) {
            new_velocity_vec[0] *= -1;
        }

        // calculate reverse matrix
        vector<vector<double> > reverse_transform(3, vector<double>(3, 0));
        reverse_transform[0][0] = y_vector[1] * z_vector[2] - y_vector[2] * z_vector[1];
        reverse_transform[0][1] = y_vector[2] * z_vector[0] - y_vector[0] * z_vector[2];
        reverse_transform[0][2] = y_vector[0] * z_vector[1] - y_vector[1] * z_vector[0];
        reverse_transform[1][0] = pos_vec[2] * z_vector[1] - pos_vec[1] * z_vector[2];
        reverse_transform[1][1] = pos_vec[0] * z_vector[2] - pos_vec[2] * z_vector[0];
        reverse_transform[1][2] = pos_vec[1] * z_vector[0] - pos_vec[0] * z_vector[1];
        reverse_transform[2][0] = pos_vec[1] * y_vector[2] - pos_vec[2] * y_vector[1];
        reverse_transform[2][1] = pos_vec[2] * y_vector[0] - pos_vec[0] * y_vector[2];
        reverse_transform[2][2] = pos_vec[0] * y_vector[1] - pos_vec[1] * y_vector[0];

        vector<double> r_velocity_vec(3, 0);
        // reverse transform to the velocity vector
        for (int i = 0; i < r_velocity_vec.size(); i++) {
            r_velocity_vec[i] = new_velocity_vec[0] * reverse_transform[i][0] + new_velocity_vec[1] * reverse_transform[i][1] + 
                new_velocity_vec[2] * reverse_transform[i][2];
        }
        
        // scale the new velocity vector to the same of original velocity
        double scale_factor(0);
        scale_factor = sqrt((pow(velocity_vec[0], 2) + pow(velocity_vec[1], 2) + pow(velocity_vec[2], 2)) / 
            (pow(r_velocity_vec[0], 2) + pow(r_velocity_vec[1], 2) + pow(r_velocity_vec[2], 2)) );
        particle.v_x = r_velocity_vec[0] * scale_factor;
        particle.v_y = r_velocity_vec[1] * scale_factor;
        particle.v_z = r_velocity_vec[2] * scale_factor;
    }
}


void Ensemble::recenter() {
    double pos_x_center(0), pos_y_center(0), pos_z_center(0);
    #pragma omp parallel for
    for (int i = 0; i < particle_number; ++i) {
        pos_x_center += ensemble[i].pos_x;
        pos_y_center += ensemble[i].pos_y;
        pos_z_center += ensemble[i].pos_z;
    }
    pos_x_center /= particle_number;
    pos_y_center /= particle_number;
    pos_z_center /= particle_number;

    // box boundary
    #pragma omp parallel for
    for (int i = 0; i < particle_number; ++i) {
        ensemble[i].pos_x += (BOX / 2 - pos_x_center);
        ensemble[i].pos_y += (BOX / 2 - pos_y_center);
        ensemble[i].pos_z += (BOX / 2 - pos_z_center);
    }
}


void Ensemble::recenter_sphere() {
    double pos_x_center(0), pos_y_center(0), pos_z_center(0);
    #pragma omp parallel for
    for (int i = 0; i < particle_number; ++i) {
        pos_x_center += ensemble[i].pos_x;
        pos_y_center += ensemble[i].pos_y;
        pos_z_center += ensemble[i].pos_z;
    }
    pos_x_center /= particle_number;
    pos_y_center /= particle_number;
    pos_z_center /= particle_number;

    // sphere boundary
    #pragma omp parallel for
    for (int i = 0; i < particle_number; ++i) {
        ensemble[i].pos_x -= pos_x_center;
        ensemble[i].pos_y -= pos_y_center;
        ensemble[i].pos_z -= pos_z_center;
    }
}


void Ensemble::rescale_temperature(double targetTemp) {
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
    double fs = sqrt(targetTemp / TEMP);
    for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
        particle->v_x = (particle->v_x - sumv_x) * fs;
        particle->v_y = (particle->v_y - sumv_y) * fs;
        particle->v_z = (particle->v_z - sumv_z) * fs;
    }
}


void Ensemble::Andersen_thermostat(double target_TEMP, double collision_frequency) {
    double sigma = sqrt(target_TEMP);
    default_random_engine random_generator;
    normal_distribution<double> gauss(0, sigma);
    uniform_real_distribution<double> ranf(0.0, 1.0);
    #pragma omp parallel for
    for (int i = 0; i < particle_number; ++i) {
        if (ranf(random_generator) < collision_frequency) {
            double scale_factor = gauss(random_generator) / calc_velocity(ensemble[i]);
            ensemble[i].v_x *= scale_factor;
            ensemble[i].v_y *= scale_factor;
            ensemble[i].v_z *= scale_factor;
        }
    }
}

double Ensemble::temperature_decreasing_curve(unsigned long i) {
    return pow(EQUILIBRATION_ITERATION - i, 2) * INIT_TEMP / EQUILIBRATION_ITERATION_2;
}


void Ensemble::energy_output(unsigned long i, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << ensemble_potential << "    " << ensemble_kinetic << "    " 
        << ensemble_potential + ensemble_kinetic << endl;
}


void Ensemble::particle_movement_output(unsigned long i, Particle& particle, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << particle.pos_x << "    " << particle.pos_y << "    " << particle.pos_z << "    " 
        << particle.v_x << "    " << particle.v_y << "    " << particle.v_z << "    " << particle.a_x_A << "    " 
        << particle.a_y_A << "    " << particle.a_z_A << "    " << particle.potential_value << "    " << particle.kinetic_value 
        << "    " << particle.potential_value + particle.kinetic_value << endl;
}


void Ensemble::temperature_output(unsigned long i, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << TEMP << endl;
}


void Ensemble::coordinates_output(ofstream& fout) {
    fout << "data_structure_1" << endl;
    fout << "_cell_length_a " << BOX << "\n_cell_length_b " << BOX << "\n_cell_length_c " << BOX << endl;
    fout << "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90" << endl;
    fout << "_cell_volume " << pow(BOX, 3) << endl;
    fout << endl;
    fout << "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z" << endl;
    for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
        fout << "O " << (it->pos_x / BOX) << " " << (it->pos_y / BOX) << " " << (it->pos_z / BOX) << endl;
    }
}



void Ensemble::inter_distance_output(unsigned long i, ofstream& fout) {
    double total_distance(0);
    #pragma omp parallel for
    for (int i = 0; i < particle_number - 1; ++i) {
        for (int j = i + 1; j < particle_number; ++j) {
            total_distance += distance(ensemble[i], ensemble[j]);
        }
    }

    fout << i * TIME_INTERVAL << "    " << total_distance << endl;
}


void Ensemble::iteration() {
    unsigned long i = 0;
    // initialize summary of velocity
    double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0), sumv2(0.0);
    while (i <= ITERATION) {
        // Initialize ensemble energy
        ensemble_kinetic = 0;
        ensemble_potential = 0;

        double temp = temperature_decreasing_curve(i);

        // calculate acceleration of step B
        #pragma omp parallel for
        for (int i = 0; i < particle_number - 1; ++i) {
            for (auto j = i + 1; j < particle_number; ++j) {
                calc_acceleration(ensemble[i], ensemble[j]);
            }
        }
        // calculate velocity of step B
        #pragma omp parallel for
        for (int i = 0; i < particle_number; ++i) {
            ensemble_potential += ensemble[i].potential_value;
            // get the instaneous temperature
            sumv_x += ensemble[i].v_x;
            sumv_y += ensemble[i].v_y;
            sumv_z += ensemble[i].v_z;
            sumv2 += ensemble[i].v_x * ensemble[i].v_x + ensemble[i].v_y * ensemble[i].v_y + ensemble[i].v_z * ensemble[i].v_z;
            // execute x and v propagation
            ensemble[i].movement();
            ensemble[i].velocity();
            if (temp > 1e-3) {
                box_boundary(ensemble[i]);
                // sphere_boundary(ensemble[i]);
            }
            // record a_A and initialize a_B
            ensemble[i].a_x_A = ensemble[i].a_x_B;
            ensemble[i].a_y_A = ensemble[i].a_y_B;
            ensemble[i].a_z_A = ensemble[i].a_z_B;
            ensemble[i].a_x_B = 0;
            ensemble[i].a_y_B = 0;
            ensemble[i].a_z_B = 0;
            ensemble[i].potential_value = 0;
            ensemble[i].kinetic_value = 0;
        }
        ensemble_kinetic = 0.5 * sumv2;
        TEMP = sumv2 / (3 * particle_number);
        sumv_x = 0;
        sumv_y = 0;
        sumv_z = 0;
        sumv2 = 0;

        if (i % SAMPLE_RATE == 0) {
            ITERATION_PERCENTAGE = ((float) i / (float) ITERATION) * 100;
            cout << "\r[MD LOG] " << get_current_time() << "\t" << ITERATION_PERCENTAGE << "\% completed " << flush;
            particle_movement_output(i, ensemble[1], particle_out);
            energy_output(i, ensemble_out);
            temperature_output(i, temperature_out);
            inter_distance_output(i, inter_distance_out);
        }

        // decrease the temperature following a parabolic curve
        if (i < EQUILIBRATION_ITERATION) {
            if (temp > 1) {
                Andersen_thermostat(temp, 0.5);
            } else {
                rescale_temperature(temp);
            }
        }

        // if (temp < 0.5 && temp > 0.49) {
        //     recenter();
        // }

        ++i;
    }
    cout << endl;   // output a new line for the progess log
    recenter();
    coordinates_output(coordinates_out);
    
    // initialize the conjugated gradient optimization
    CG conjugated_gradient_optimization(outputPath, ensemble);
    conjugated_gradient_optimization.conjugated_gradient_minimization(BOX);
}


#endif
