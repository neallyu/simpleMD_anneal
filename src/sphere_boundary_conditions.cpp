#include <iostream>
#include <cmath>
#include <vector>
#include <iterator>
#include <algorithm>

using namespace std;

int main() {
    vector<float> pos_vec({1, 0 ,0});
    vector<float> velocity_vec({-1, 1, 1});

    float a = -(pow(pos_vec[0], 2) + pow(pos_vec[1], 2) + pow(pos_vec[2], 2)) / 
        (pos_vec[0] * velocity_vec[0] + pos_vec[1] * velocity_vec[1] + pos_vec[2] * velocity_vec[2]);

    vector<float> y_vector(pos_vec);
    y_vector[0] *= a; y_vector[0] += velocity_vec[0];
    y_vector[1] *= a; y_vector[1] += velocity_vec[1];
    y_vector[2] *= a; y_vector[2] += velocity_vec[2];

    vector<float> z_vector(3, 0);
    z_vector[0] = pos_vec[1] * y_vector[2] - pos_vec[2] * y_vector[1];
    z_vector[1] = pos_vec[2] * y_vector[0] - pos_vec[0] * y_vector[2];
    z_vector[2] = pos_vec[0] * y_vector[1] - pos_vec[1] * y_vector[0];

    // linear transformation of velocity vector from original cooridinates to new cooridnates
    vector<float> new_velocity_vec(3, 0);
    new_velocity_vec[0] = pos_vec[0] * velocity_vec[0] + y_vector[0] * velocity_vec[1] + z_vector[0] * velocity_vec[2];
    new_velocity_vec[1] = pos_vec[1] * velocity_vec[0] + y_vector[1] * velocity_vec[1] + z_vector[1] * velocity_vec[2];
    new_velocity_vec[2] = pos_vec[2] * velocity_vec[0] + y_vector[2] * velocity_vec[1] + z_vector[2] * velocity_vec[2];

    // if the x-portion of the vector larger than zero, mirror operation of transformed velocity vector correpsonding x axis
    if (new_velocity_vec[0] > 0) {
        new_velocity_vec[0] *= -1;
    }

    // calculate reverse matrix
    vector<vector<float> > reverse_transform(3, vector<float>(3, 0));
    reverse_transform[0][0] = y_vector[1] * z_vector[2] - y_vector[2] * z_vector[1];
    reverse_transform[0][1] = y_vector[2] * z_vector[0] - y_vector[0] * z_vector[2];
    reverse_transform[0][2] = y_vector[0] * z_vector[1] - y_vector[1] * z_vector[0];
    reverse_transform[1][0] = pos_vec[2] * z_vector[1] - pos_vec[1] * z_vector[2];
    reverse_transform[1][1] = pos_vec[0] * z_vector[2] - pos_vec[2] * z_vector[0];
    reverse_transform[1][2] = pos_vec[1] * z_vector[0] - pos_vec[0] * z_vector[1];
    reverse_transform[2][0] = pos_vec[1] * y_vector[2] - pos_vec[2] * y_vector[1];
    reverse_transform[2][1] = pos_vec[2] * y_vector[0] - pos_vec[0] * y_vector[2];
    reverse_transform[2][2] = pos_vec[0] * y_vector[1] - pos_vec[1] * y_vector[0];

    vector<float> r_velocity_vec(3, 0);
    // reverse transform to the velocity vector
    for (int i = 0; i < r_velocity_vec.size(); i++) {
        r_velocity_vec[i] = new_velocity_vec[0] * reverse_transform[i][0] + new_velocity_vec[1] * reverse_transform[i][1] + 
            new_velocity_vec[2] * reverse_transform[i][2];
    }
    
    // scale the new velocity vector to the same of original velocity
    float scale_factor(0);
    scale_factor = sqrt((pow(velocity_vec[0], 2) + pow(velocity_vec[1], 2) + pow(velocity_vec[2], 2)) / 
        (pow(r_velocity_vec[0], 2) + pow(r_velocity_vec[1], 2) + pow(r_velocity_vec[2], 2)) );
    velocity_vec[0] = r_velocity_vec[0] * scale_factor;
    velocity_vec[1] = r_velocity_vec[1] * scale_factor;
    velocity_vec[2] = r_velocity_vec[2] * scale_factor;

    copy(velocity_vec.begin(), velocity_vec.end(), ostream_iterator<float>(cout, " ")); cout << endl;
}