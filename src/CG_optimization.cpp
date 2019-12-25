#ifndef CG_H
#define CG_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <string>
// #include <omp.h>
#include "anneal/utils.h"
#include "anneal/particle.h"
#include "anneal/error_handling.h"
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

class CG {

public:
    CG(string _output_path, MatrixXd& ensemble);
    ~CG();

    inline MatrixXd calc_gradient(MatrixXd& position);

    inline double calc_potential(MatrixXd& position);

    // calculate d_ii = -g_ii + beta * d_i
    inline MatrixXd update_direction(MatrixXd& d_i, MatrixXd& g_ii, MatrixXd& g_i);

    // calculate x_ii = x_i + alpha * d_i
    inline MatrixXd line_search(MatrixXd& position, MatrixXd& g_i, MatrixXd& d_i);

    inline void conjugated_gradient_minimization(const double BOX);

    inline void trend_output(double potential_value, double gradient_squaredNorm, ofstream& fout);

    inline void coordinates_output(const double BOX, ofstream& fout);

private:
    MatrixXd position;                              // store the position of the particles
    // vector<vector<double> > position;               // store the position of the particles
    double c1;                                      // constant in Wolfe line search
    double c2;                                      // constant in Wolfe line search 0 < c1 < c2 < 1
    ofstream trend_out;                             // output file stream of potential gradient norm
    ofstream coordinates_out_CG;                    // output file stream of coordinates after CG optimization
};


CG::CG(string _output_path, MatrixXd& ensemble): 
    position(ensemble),
    c1(1e-4),
    c2(0.8),
    coordinates_out_CG(_output_path + "/coordinates_after_CG.cif"),
    trend_out(_output_path + "/trend.csv")
    { }


CG::~CG() {
    coordinates_out_CG.close();
    trend_out.close();
}


MatrixXd CG::calc_gradient(MatrixXd& position) {
    MatrixXd gradient(position.rows(), 3);
    gradient.fill(0);
    // #pragma omp parallel for
    for (int i = 0; i < position.rows() - 1; ++i) {
        for (int j = i + 1; j < position.rows(); ++j) {
            double dx = position(i, 0) - position(j, 0);
            double dy = position(i, 1) - position(j, 1);
            double dz = position(i, 2) - position(j, 2);
            double r2 = dx * dx + dy * dy + dz * dz;
            double r2i = 1 / r2;
            double r6i = pow(r2i, 3);
            double g_x = -48.0 * r6i * r2i * (r6i - 0.5) * dx;
            double g_y = -48.0 * r6i * r2i * (r6i - 0.5) * dy;
            double g_z = -48.0 * r6i * r2i * (r6i - 0.5) * dz;
            gradient(i, 0) += g_x;
            gradient(i, 1) += g_y;
            gradient(i, 2) += g_z;
            gradient(j, 0) -= g_x;
            gradient(j, 1) -= g_y;
            gradient(j, 2) -= g_z;
        }
    }
    return gradient;
}


double CG::calc_potential(MatrixXd &position) {
    double potential_value(0);
    // #pragma omp parallel for
    for (int i = 0; i < position.rows() - 1; ++i) {
        for (int j = i + 1; j < position.rows(); ++j) {
            double dx = position(i, 0) - position(j, 0);
            double dy = position(i, 1) - position(j, 1);
            double dz = position(i, 2) - position(j, 2);
            double r2 = dx * dx + dy * dy + dz * dz;
            double r2i = 1 / r2;
            double r6i = pow(r2i, 3);
            potential_value += 4.0 * r6i * (r6i - 1);
        }
    }
    return potential_value;
}


MatrixXd CG::update_direction(MatrixXd& d_i, MatrixXd& g_ii, MatrixXd& g_i) {
    double beta = g_ii.squaredNorm() / g_i.squaredNorm();
    MatrixXd d_ii = d_i * beta - g_ii;
    return d_ii;
}


MatrixXd CG::line_search(MatrixXd& position, MatrixXd& g_i, MatrixXd& d_i) {
    // calculate the constant in the line search of this time
    double alpha_step_size(1e-7);
    double alpha(1e-3);
    double potential = calc_potential(position);
    // the sum of diagonal elements, which is x1 * d_x1 + y1 * d_y1 + z1 * d_z1 + ...
    double constant = (g_i.transpose() * d_i).trace();
    double constant_abs = fabs(constant);

    // the mutable part in the line search of this time
    MatrixXd position_search = position + alpha * d_i;
    double position_search_potential = calc_potential(position_search);
    double position_search_variable_abs = fabs((calc_gradient(position_search).transpose() * d_i).trace());

    while (position_search_potential > potential + alpha * c1 * constant || position_search_variable_abs > c2 * constant_abs) {
        alpha -= alpha_step_size;
        if (alpha <= 0) {
            alpha = alpha_step_size;
            position_search = position + alpha * d_i;
            break;
        }
        position_search = position + alpha * d_i;
        position_search_potential = calc_potential(position_search);
        position_search_variable_abs = fabs((calc_gradient(position_search).transpose() * d_i).trace());
    }
    return position_search;
}


void CG::conjugated_gradient_minimization(const double BOX) {
    cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization started..." << endl;
    // initialize
    double precision_2(1e-3);
    MatrixXd position_i(position);
    MatrixXd g_i = calc_gradient(position);
    MatrixXd d_i = -1 * g_i;
    MatrixXd position_ii(position_i);
    MatrixXd d_ii(d_i);
    MatrixXd g_ii(g_i);
    double potential_value = calc_potential(position_ii);
    double gradient_squaredNorm = g_ii.squaredNorm();

    cout << "potential value: " << potential_value << endl;
    cout << "gradient squared norm: " << gradient_squaredNorm << endl;
    cout << "precision2: " << precision_2 << endl;
    trend_output(potential_value, gradient_squaredNorm, trend_out);

    // minimal condition is the norm of gradient of fx closes to zero
    while (gradient_squaredNorm > precision_2) {
        position_ii = line_search(position_i, g_i, d_i);
        position_i = position_ii;
        g_ii = calc_gradient(position_ii);
        d_ii = update_direction(d_i, g_ii, g_i);
        d_i = d_ii;
        g_i = g_ii;
        gradient_squaredNorm = g_ii.squaredNorm();
        potential_value = calc_potential(position_ii);
        trend_output(potential_value, gradient_squaredNorm, trend_out);
    }

    //record final position
    position = position_ii;
    coordinates_output(BOX, coordinates_out_CG);
    cout << "[MD LOG] " << get_current_time() << "\tConjugated gradient optimization finished." << endl;
}


void CG::trend_output(double potential_value, double gradient_squaredNorm, ofstream& fout) {
    fout << potential_value << "    " << gradient_squaredNorm << endl;
}


void CG::coordinates_output(const double BOX, ofstream& fout) {
    fout << "data_structure_1" << endl;
    fout << "_cell_length_a " << BOX << "\n_cell_length_b " << BOX << "\n_cell_length_c " << BOX << endl;
    fout << "_cell_angle_alpha 90\n_cell_angle_beta 90\n_cell_angle_gamma 90" << endl;
    fout << "_cell_volume " << pow(BOX, 3) << endl;
    fout << endl;
    fout << "loop_\n_atom_site_label\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z" << endl;
    for (int i = 0; i < position.rows(); i++) {
        fout << "O " << position(i, 0) / BOX << " " << position(i, 0) / BOX << " " << position(i, 2) / BOX << endl;
    }
}

#endif


int main(int argc, char* argv[]) {
    int BOX = 0;
    int particleNumber = 0;
    // omp_set_num_threads(32);
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

    MatrixXd ensemble(particleNumber, 3);
    for (int i = 0; i < particleNumber; i++) {
        string coordinateX, coordinateY, coordinateZ;
        input >> coordinateX >> coordinateY >> coordinateZ;
        ensemble(i, 0) = stod(coordinateX);
        ensemble(i, 1) = stod(coordinateY);
        ensemble(i, 2) = stod(coordinateZ);
    }
    
    CG cg_optimization(output_path, ensemble);
    cg_optimization.conjugated_gradient_minimization(BOX);
}
