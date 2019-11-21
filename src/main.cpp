#include <fstream>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include "error_handling.hpp"
#include "ensemble_reduced_unit.hpp"
// #include "ensemble_reduced_unit_for_msd.hpp"
#include "utils.hpp"

using namespace std;

void read_input_reduced_unit(int argc, char *argv[]){
    unsigned particle_number;
    double init_temp;
    double set_temp;
    double time_interval;
    double equilibration_time;
    double total_time;
    double rho;

    try {
        ifstream input;
        input.open(argv[1]);
        if (!input) {
            throw ReadingFile_Open();
        }
        cout << "[MD LOG] " << get_current_time() << "\tReading the input file \"" << argv[1] << "\"..." << endl;
        int i = 0;
        string parameter;
        while (input >> parameter) {
            switch (i) {
                case 1:
                    particle_number = stoi(parameter);
                    break;
                case 3:
                    init_temp = stod(parameter);
                    break;
                case 5:
                    set_temp = stod(parameter);
                    break;
                case 7:
                    time_interval = stod(parameter);
                    break;
                case 9:
                    equilibration_time = stod(parameter);
                    break;
                case 11:
                    total_time = stod(parameter);
                    break;
                case 13:
                    rho = stod(parameter);
                    break;
            }
            ++i;
        }
        if (i != 14) {
            throw ReadingFile_Other();
        }
        if (!mkdir(argv[2], S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) {
            throw CreatingOutputPath();
        }

    } catch (ReadingFile_Open e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
    } catch (ReadingFile_Other e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
    } catch (CreatingOutputPath e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in: " << e.what() << endl;
    }

    cout << "[MD LOG] " << get_current_time() << "\tParameters is successfully inputed" << endl;
    cout << "[MD LOG] " << get_current_time() << "\tParticle number: " << particle_number << endl;
    cout << "[MD LOG] " << get_current_time() << "\tInitial temperature: " << init_temp << endl;
    cout << "[MD LOG] " << get_current_time() << "\tEquilibration temperature: " << set_temp << endl;
    cout << "[MD LOG] " << get_current_time() << "\tEquilibration time: " << equilibration_time << endl;
    cout << "[MD LOG] " << get_current_time() << "\tTotal time: " << total_time << endl;
    cout << "[MD LOG] " << get_current_time() << "\tDensity: " << rho << endl;

    cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    
    // Initialize the ensemble
    Ensemble ensemble1(particle_number, init_temp, set_temp, time_interval, equilibration_time, total_time, rho, argv[2]);

    cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
    ensemble1.iteration();

    cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
}


// void read_input_real_unit(int argc, char* argv[]) {
//     string molecule;
//     unsigned particle_number;
//     double sigma;
//     double epsilon;
//     double MASS;
//     double init_temp;
//     double set_temp;
//     double time_interval;
//     double equilibration_time;
//     double total_time;
//     double box;

//     try {
//         ifstream input;
//         input.open(argv[1]);
//         if (!input) {
//             throw ReadingFile_Open();
//         }
//         cout << "[MD LOG] " << get_current_time() << "\tReading the input file \"" << argv[1] << "\"..." << endl;
//         int i = 0;
//         string parameter;
//         while (input >> parameter) {
//             switch (i) {
//                 case 1:
//                     molecule = parameter;
//                     break;
//                 case 3:
//                     particle_number = stoi(parameter);
//                     break;
//                 case 5:
//                     sigma = stod(parameter);
//                     break;
//                 case 7:
//                     epsilon = stod(parameter);
//                     break;
//                 case 9:
//                     MASS = stod(parameter);
//                     break;            
//                 case 11:
//                     init_temp = stod(parameter);
//                     break;
//                 case 13:
//                     set_temp = stod(parameter);
//                     break;
//                 case 15:
//                     time_interval = stod(parameter);
//                     break;
//                 case 17:
//                     equilibration_time = stod(parameter);
//                     break;
//                 case 19:
//                     total_time = stod(parameter);
//                     break;
//                 case 21:
//                     box = stod(parameter);
//                     break;
//             }
//             ++i;
//         }
//         if (i != 22) {
//             throw ReadingFile_Other();
//         }
//         if (!mkdir(argv[2], S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) {
//             throw CreatingOutputPath();
//         }

//     } catch (ReadingFile_Open e) {
//         cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
//     } catch (ReadingFile_Other e) {
//         cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
//     } catch (CreatingOutputPath e) {
//         cerr << "[MD ERR] " << get_current_time() << "\tError in: " << e.what() << endl;
//     }

//     cout << "[MD LOG] " << get_current_time() << "\tParameters is successfully inputed" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tmolecule: " << molecule << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tparticle number: " << particle_number << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tsigma: " << sigma << " Angstrom" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tepsilon: " << epsilon << " kJ/mol"<< endl;
//     cout << "[MD LOG] " << get_current_time() << "\tmass: " << MASS << " g/mol" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tinitial temperature: " << init_temp << " K" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tequilibration temperature: " << set_temp << " K" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\ttime interval: " << time_interval << " fs" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tequilibration time: " << equilibration_time << " ns" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\ttotal time: " << total_time << " ns" << endl;
//     cout << "[MD LOG] " << get_current_time() << "\tbox size: " << box << " Angstrom" << endl;

//     cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    
//     // Initialize the ensemble
//     Ensemble ensemble1(particle_number, sigma, epsilon, MASS, init_temp, set_temp, time_interval, equilibration_time, total_time, box, argv[2]);

//     cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
//     ensemble1.iteration();

//     cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
// }


int main(int argc, char *argv[]) {
    read_input_reduced_unit(argc, argv);
}
