#ifndef CG_H
#define CG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include <omp.h>
#include "utils.hpp"
#include "particle.hpp"
#include "error_handling.hpp"

using namespace std;

class CG {

public:
    CG(string _output_path, vector<vector<double> >& ensemble);
    ~CG();

    inline vector<vector<double> > calc_potential_gradient(vector<vector<double> >& position);

    inline double calc_potential_value(vector<vector<double> > &position);

    inline double calc_potential_gradient_norm_2(vector<vector<double> >& gradient);

    inline void calc_direction(vector<vector<double> >& d_i, 
        vector<vector<double> >& gradient_ii, vector<vector<double> >& gradient_i);

    inline vector<vector<double> > matrix_addition_x(vector<vector<double> >& position, vector<vector<double> >& d_i, double& alpha);
    inline vector<vector<double> > matrix_addition_d(vector<vector<double> >& g_ii, vector<vector<double> >& d_i, double& beta);

    inline vector<vector<double> > line_search(vector<vector<double> >& position, vector<vector<double> >& g_i, vector<vector<double> >& d_i);

    inline void conjugated_gradient_minimization(const double BOX);

    inline void potential_value_output(double value, ofstream& fout);

    inline void coordinates_output(const double BOX, ofstream& fout);

private:
    vector<vector<double> > position;               // store the position of the particles
    double c;                                       // constant in Wolfe line search
    ofstream potential_value_out;                   // output file stream of potential gradient norm
    ofstream coordinates_out_CG;                    // output file stream of coordinates after CG optimization
};


CG::CG(string _output_path, vector<vector<double> >& ensemble): 
    position(ensemble.size(), vector<double>(3, 0)),
    c(1e-4),
    coordinates_out_CG(_output_path + "/coordinates_after_CG.cif"),
    potential_value_out(_output_path + "/potential_value.csv")
    {
        // parallel
        omp_set_num_threads(32);

        for (int i = 0; i < ensemble.size(); i++) {
            position[i][0] = ensemble[i][0];
            position[i][1] = ensemble[i][1];
            position[i][2] = ensemble[i][2];
        }
}


CG::~CG() {
    coordinates_out_CG.close();
    potential_value_out.close();
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


void CG::calc_direction(vector<vector<double> >& d_i,
    vector<vector<double> >& gradient_ii, vector<vector<double> >& gradient_i) {
    
    double beta = calc_potential_gradient_norm_2(gradient_ii) / calc_potential_gradient_norm_2(gradient_i);

    #pragma omp parallel for
    for (int i = 0; i < d_i.size(); i++) {
        d_i[i][0] *= beta;
        d_i[i][0] -= gradient_ii[i][0];
        d_i[i][1] *= beta;
        d_i[i][1] -= gradient_ii[i][1];
        d_i[i][2] *= beta;
        d_i[i][2] -= gradient_ii[i][2];
    }
}


vector<vector<double> > CG::matrix_addition_x(vector<vector<double> >& position, vector<vector<double> >& d_i, double& alpha) {
    vector<vector<double> > result(position.size(), vector<double>(3, 0));
    #pragma omp parallel for
    for (int i = 0; i < position.size(); i++) {
        result[i][0] = position[i][0] + alpha * d_i[i][0];
        result[i][1] = position[i][1] + alpha * d_i[i][1];
        result[i][2] = position[i][2] + alpha * d_i[i][2];
    }

    return result;
}


vector<vector<double> > CG::matrix_addition_d(vector<vector<double> >& g_ii, vector<vector<double> >& d_i, double& beta) {
    vector<vector<double> > result(g_ii.size(), vector<double>(3, 0));
    #pragma omp parallel for
    for (int i = 0; i < g_ii.size(); i++) {
        result[i][0] = -g_ii[i][0] + beta * d_i[i][0];
        result[i][1] = -g_ii[i][1] + beta * d_i[i][1];
        result[i][2] = -g_ii[i][2] + beta * d_i[i][2];
    }

    return result;
}


vector<vector<double> > CG::line_search(vector<vector<double> >& position, vector<vector<double> >& g_i, vector<vector<double> >& d_i) {
    double alpha(1e-8);
    double constant(0);

    #pragma omp parallel for
    for (int i = 0; i < g_i.size(); i++) {
        constant += g_i[i][0] * d_i[i][0] + g_i[i][1] * d_i[i][1] + g_i[i][2] * d_i[i][2];
    }

    vector<vector<double> > position_search = matrix_addition_x(position, d_i, alpha);
    double potential_value = calc_potential_value(position);
    double potential_value_search = calc_potential_value(position_search);

    while (potential_value_search > potential_value + alpha * c * constant) {
        alpha += 1e-8;
        position_search = matrix_addition_x(position, d_i, alpha);
        potential_value_search = calc_potential_value(position_search);
    }

    return position_search;
}


void CG::conjugated_gradient_minimization(const double BOX) {
    cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization started..." << endl;

    // final precision
    double precision_2(1e-8);

    // initialize position matrix
    vector<vector<double> > position_i(position);

    // initialize g_i
    vector<vector<double> > g_i = calc_potential_gradient(position);

    // initialize d_i
    vector<vector<double> > d_i(g_i);
    #pragma omp parallel for
    for (int i = 0; i < g_i.size(); i++) {
        d_i[i][0] *= -1;
        d_i[i][1] *= -1;
        d_i[i][2] *= -1;
    }

    // intialize position_ii, d_ii and g_ii
    vector<vector<double> > position_ii(position_i);
    vector<vector<double> > d_ii(d_i);
    vector<vector<double> > g_ii(g_i);

    //initalize potential value
    double potential_value = calc_potential_value(position_ii);

    // minimal condition is the norm of gradient of fx closes to zero
    while (potential_value > precision_2) {
        potential_value_output(potential_value, potential_value_out);

        // calculate x_ii = x_i + alpha * d_i
        position_ii = line_search(position_i, g_i, d_i);

        // calculate d_ii = -g_ii + beta * d_i
        g_ii = calc_potential_gradient(position_ii);
        double beta = calc_potential_gradient_norm_2(g_ii) / calc_potential_gradient_norm_2(g_i);
        d_ii = matrix_addition_d(g_ii, d_i, beta);

        // update potential value
        potential_value = calc_potential_value(position_ii);
    }

    //record final position
    for (int i = 0; i < position_ii.size(); i++) {
        position[i][0] = position_ii[i][0];
        position[i][1] = position_ii[i][1];
        position[i][2] = position_ii[i][2];
    }

    coordinates_output(BOX, coordinates_out_CG);
    cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization finished." << endl;
}


void CG::potential_value_output(double value, ofstream& fout) {
    fout << value << endl;
}


void CG::coordinates_output(const double BOX, ofstream& fout) {
    fout << "data_structure_1" << endl;
    fout << "_cell_length_a " << BOX << "\n_cell_length_b " << BOX << "\n_cell_length_c " << BOX << endl;
    fout << "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90" << endl;
    fout << "_cell_volume " << pow(BOX, 3) << endl;
    fout << endl;
    fout << "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z" << endl;
    for (int i = 0; i < position.size(); i++) {
        fout << "O " << position[i][0] / BOX << " " << position[i][1] / BOX << " " << position[i][2] / BOX << endl;
    }
}

#endif


int main(int argc, char* argv[]) {
    int BOX = 0;
    int particleNumber = 0;
        omp_set_num_threads(32);
    ifstream input;
    string output_path;
    try {
        input.open(argv[1]);
        if (!input) {
            throw ReadingFile_Open();
        }
        cout << "[MD LOG] " << get_current_time() << "\tReading the input file \"" << argv[1] << "\"..." << endl;
        BOX = stoi(argv[2]);
        particleNumber = stoi(argv[3]);
        if (BOX <= 0 || particleNumber <= 0) {
            throw WrongInputNumber();
        }
        output_path = argv[4];
    } catch (ReadingFile_Open e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
    } catch (WrongInputNumber e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in getting BOX size: " << BOX << " or in particle number: " << particleNumber << endl;
    }

    vector<vector<double> > ensemble(particleNumber, vector<double>(3, 0));
    for (int i = 0; i < particleNumber; i++) {
        string coordinateX, coordinateY, coordinateZ;
        input >> coordinateX >> coordinateY >> coordinateZ;
        ensemble[i][0] = stod(coordinateX) * BOX;
        ensemble[i][1] = stod(coordinateY) * BOX;
        ensemble[i][2] = stod(coordinateZ) * BOX;
    }
    
    CG cg_optimization(output_path, ensemble);
    cg_optimization.conjugated_gradient_minimization(BOX);
}