#ifndef CG_H
#define CG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>
#include "utils.hpp"

using namespace std;

class CG {
public:
    CG(string _output_path);

    inline vector<vector<double> > calc_potential_gradient(vector<vector<double> >& position);

    inline double calc_potential_value(vector<vector<double> > &position);

    inline double calc_potential_gradient_norm_2(vector<vector<double> >& gradient);

    inline double line_search(vector<vector<double> >& position);

    // inline vector<double> calc_gamma_i();

    inline void conjugated_gradient_minimization();

    inline void potential_gradient_norm_output(int count, double norm, ofstream& fout);
private:
    ofstream potential_gradient_norm_out;           // output file stream of potential gradient norm
    ofstream coordinates_out_CG;                    // output file stream of coordinates after CG optimization
};


CG::CG(string _output_path): 
    coordinates_out_CG(_output_path + "/coordinates_after_CG.cif"),
    potential_gradient_norm_out(_output_path + "/potential_gradient_norm.csv") 
    {
        // parallel
        omp_set_num_threads(32);
}


CG::~CG() {
    coordinates_out_CG.close();
    potential_gradient_norm_out.close();
}


vector<vector<double> > CG::calc_potential_gradient(vector<vector<double> >& position) {

    vector<vector<double> > gradient(position.size(), vector<double>(3));

    #pragma omp parallel for
    for (int i = 0; i < position.size() - 1; ++i) {
        for (int j = i + 1; j < position.size(); ++j) {
            double dx = position[i][0] - position[j][0];
            double dy = position[i][1] - position[j][1];
            double dz = position[i][2] - position[j][2];

            double r2 = dx * dx + dy * dy + dz * dz;
            double r2i = 1 / r2;
            double r6i = pow(r2i, 3);
            double g_x = -48.0 * r6i * r2i * (r6i - 0.5) * dx;
            double g_y = -48.0 * r6i * r2i * (r6i - 0.5) * dy;
            double g_z = -48.0 * r6i * r2i * (r6i - 0.5) * dz;
            gradient[i][0] += g_x;
            gradient[i][1] += g_y;
            gradient[i][2] += g_z;
            gradient[j][0] -= g_x;
            gradient[j][1] -= g_y;
            gradient[j][2] -= g_z;
        }
    }

    return gradient;
}


double CG::calc_potential_value(vector<vector<double> > &position) {
    double potential_value(0);

    #pragma omp parallel for
    for (int i = 0; i < position.size() - 1; ++i) {
        for (int j = i + 1; j < position.size(); ++j) {
            double dx = position[i][0] - position[j][0];
            double dy = position[i][1] - position[j][1];
            double dz = position[i][2] - position[j][2];

            double r2 = dx * dx + dy * dy + dz * dz;
            double r2i = 1 / r2;
            double r6i = pow(r2i, 3);
            potential_value += 4.0 * r6i * (r6i - 1);
        }
    }

    return potential_value;
}


double CG::calc_potential_gradient_norm_2(vector<vector<double> >& gradient) {
    double norm_2(0);

    #pragma omp parallel for
    for (int i = 0; i < gradient.size(); ++i) {
        norm_2 += pow(gradient[i][0], 2) + pow(gradient[i][1], 2) + pow(gradient[i][2], 2);
    }
    return norm_2;
}


// vector<double> Ensemble::calc_gamma_i() {
//     vector<double> gamma_i(2, 0);
//     #pragma omp parallel for
//     for (int i = 0; i < particle_number; ++i) {
//         gamma_i[0] += pow(ensemble[i].a_x_A, 2) + pow(ensemble[i].a_y_A, 2) + pow(ensemble[i].a_z_A, 2);
//         gamma_i[1] += (ensemble[i].a_x_B - ensemble[i].a_x_A) * ensemble[i].a_x_B +
//             (ensemble[i].a_y_B - ensemble[i].a_y_A) * ensemble[i].a_y_B +
//             (ensemble[i].a_z_B - ensemble[i].a_z_A) * ensemble[i].a_z_B;
//     }
//     return gamma_i;
// }


double CG::line_search(vector<vector<double> >& position) {
    vector<vector<double> > position_search(position);
    double alpha(1e-4);
    


}


// void CG::conjugated_gradient_minimization() {
//     cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization started..." << endl;

//     // declare ensemble for conjugated gradient optimization
//     // vector<Particle_CG> ensemble_CG;
//     double step_size(1e-60), precision_2(1e-2);

//     vector<vector<int> > pos_i(particle_number, vector<int>(3));

//     // initialize pos_i 2D vector
//     for (int i = 0; i < particle_number; ++i) {
//         pos_i[i][0] = ensemble[i].pos_x;
//         pos_i[i][1] = ensemble[i].pos_y;
//         pos_i[i][2] = ensemble[i].pos_z;
//     }
    


//     // initialize d_i
//     #pragma omp parallel for
//     for (int i = 0; i < particle_number - 1; ++i) {
//         for (auto j = i + 1; j < particle_number; ++j) {
//             calc_acceleration(ensemble[i], ensemble[j]);
//         }
//     }
//     #pragma omp parallel for
//     for (int i = 0; i < particle_number; ++i) {
//         ensemble[i].d_x_A = ensemble[i].a_x_B;
//         ensemble[i].d_y_A = ensemble[i].a_y_B;
//         ensemble[i].d_z_A = ensemble[i].a_z_B;
//     }

//     vector<double> g = calc_potential_gradient_norm_2();
//     cout << "gradient_norm_2: " << g[0] << endl;
//     int count(0);

//     // minimal condition is the norm of gradient of fx closes to zero
//     while (g[0] > precision_2) {

//         potential_gradient_norm_output(count, g[0], potential_gradient_norm_out);

//         #pragma omp parallel for
//         for (int i = 0; i < particle_number; ++i) {
//             // x propagation: x_(i+1) = x_i + step_size * d_i
//             ensemble[i].pos_x += step_size * ensemble[i].d_x_A;
//             ensemble[i].pos_y += step_size * ensemble[i].d_y_A;
//             ensemble[i].pos_z += step_size * ensemble[i].d_z_A;

//             // record a_A and initialize a_B
//             ensemble[i].a_x_A = ensemble[i].a_x_B;
//             ensemble[i].a_y_A = ensemble[i].a_y_B;
//             ensemble[i].a_z_A = ensemble[i].a_z_B;
//             ensemble[i].a_x_B = 0;
//             ensemble[i].a_y_B = 0;
//             ensemble[i].a_z_B = 0;
//         }

//         // calculate g_(i+1) 
//         #pragma omp parallel for
//         for (int i = 0; i < particle_number - 1; ++i) {
//             for (auto j = i + 1; j < particle_number; ++j) {
//                 calc_acceleration(ensemble[i], ensemble[j]);
//             }
//         }

//         // g = calc_potential_gradient_norm_2();
//         // double beta_i = g[1] / g[0];

//         g = calc_gamma_i();
//         double beta_i = g[1] / g[0];
        
//         // d propagation: d_(i+1) = a_(i+1) + norm_2(a_(i+1)) / norm_2(a_i) * d_i
//         #pragma omp parallel for
//         for (int i = 0; i < particle_number; ++i) {
//             ensemble[i].d_x_A = ensemble[i].d_x_B;
//             ensemble[i].d_y_A = ensemble[i].d_y_B;
//             ensemble[i].d_z_A = ensemble[i].d_z_B;
//             ensemble[i].d_x_B = ensemble[i].a_x_B + beta_i * ensemble[i].d_x_A;
//             ensemble[i].d_y_B = ensemble[i].a_y_B + beta_i * ensemble[i].d_y_A;
//             ensemble[i].d_z_B = ensemble[i].a_z_B + beta_i * ensemble[i].d_z_A;
//         }
//         count++;
//     }

//     coordinates_output(coordinates_out_CG);
//     cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization finished." << endl;
// }




// void Ensemble::potential_gradient_norm_output(int count, double norm, ofstream& fout) {
//     fout << count << "    " << norm << "    " << endl;
// }


#endif