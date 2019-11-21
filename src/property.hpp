#ifndef PROPERTY_H
#define PROPERTY_H

#include "particle.hpp"
#include "utils.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class Property {

public:
    Property(int _particle_number, double _time_interval, string _output_path):
        v_index(0), autocorr(0), v_avg(0), TIME_INTERVAL(_time_interval), particle_number(_particle_number),
        velocity_autocorr_out(_output_path + "/velocity_autocorr.csv") { }

    void initalize(vector<Particle>& ensemble) {
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Start_status.push_back(*it);
        }
    }

    double calc_mean_square_particle_displacement(vector<Particle>& ensemble) {
        double MSD(0);
        for (int i = 0; i < ensemble.size(); ++i) {
            MSD += (pow((ensemble[i].pos_x - Start_status[i].pos_x), 2) 
                + pow((ensemble[i].pos_y - Start_status[i].pos_y), 2)
                + pow((ensemble[i].pos_z - Start_status[i].pos_z), 2));
        }
        return (MSD / (double) particle_number);
    }

    // record v_x of each particle at the moment
    void sample_velocity_autocorrelation(vector<Particle> &ensemble) {
        Velocity.resize(v_index + 1);
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Velocity[v_index].push_back(*it);
            v_avg += it->v_x;
        }
        ++v_index;
    }


    void calc_velocity_autocorrelation() {
        cout << "\r[MD LOG] " << get_current_time() << "\t" << "Calculating velocity autocorrelation function" << endl;
        double ITERATION_PERCENTAGE(0);
        double V_VARIANCE(0);
        v_avg /= (v_index * particle_number);
        // iteration counter
        unsigned long I = 0;
        unsigned long TOTAL = v_index * (v_index + 1) * particle_number;
        // #pragma omp parallel for
        for (int t_frame = 0; t_frame < v_index; ++t_frame) {
            for (int particle_index = 0; particle_index < Velocity[0].size(); ++ particle_index) {
                for (int n = 0; n < v_index - t_frame; ++n) {
                    autocorr += ((Velocity[n][particle_index].v_x - v_avg) * (Velocity[n + t_frame][particle_index].v_x - v_avg));
                    ITERATION_PERCENTAGE = ((double) I / (double) TOTAL) * 100;
                    I++;
                    cout << "\r[MD LOG] " << get_current_time() << "\t" << ITERATION_PERCENTAGE << "\% completed       " << flush;
                }
                V_VARIANCE += pow(Velocity[t_frame][particle_index].v_x - v_avg, 2);
            }
            velocity_autocorr.push_back(autocorr);
            autocorr = 0;

        }
        // cout << endl;
        for (auto it = velocity_autocorr.begin(); it != velocity_autocorr.end(); ++it) {
            *it /= V_VARIANCE;
        }
    }


    void velocity_autocorr_output() {
        for (int i = 0; i < velocity_autocorr.size(); ++i){
            velocity_autocorr_out << i * TIME_INTERVAL << "    " << velocity_autocorr[i] << endl;
        }
        velocity_autocorr_out.close();
    }


private:
    vector<Particle> Start_status;
    vector<vector <Particle> > Velocity;
    int v_index;
    int particle_number;
    double v_avg;
    double autocorr;
    vector<double> velocity_autocorr;
    double TIME_INTERVAL;
    ofstream velocity_autocorr_out;
};

#endif