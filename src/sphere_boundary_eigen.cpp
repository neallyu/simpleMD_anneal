#include <iostream>
#include "Eigen/Dense"

using namespace std;
using namespace Eigen;

int main() {
    Vector3f x_vec(0, 0, 1);
    Vector3f v_vec(1, 1, 1);
    cout << "position:\n" << x_vec << endl;
    cout << "velocity:\n" << v_vec << endl;
    float a = -x_vec.squaredNorm() / (x_vec.dot(v_vec));

    Vector3f y_vec(x_vec);
    y_vec *= a;
    y_vec += v_vec;

    Vector3f z_vec(0, 0, 0);
    z_vec = x_vec.cross(y_vec);

    Matrix3f transform_matrix_transpose, transform_matrix;
    transform_matrix_transpose << x_vec, y_vec, z_vec;
    transform_matrix = transform_matrix_transpose.transpose();
    cout << "transform_matrix:\n" << transform_matrix << endl;

    Vector3f transformed_v_vec(0, 0, 0);
    transformed_v_vec = transform_matrix * v_vec;
    cout << "tranformed_v_vec:\n" << transformed_v_vec << endl;

    if (transformed_v_vec[0] > 0) {
        transformed_v_vec[0] *= -1;
    }

    Matrix3f reverse_transform_matrix = transform_matrix.inverse();
    Vector3f reverse_transformed_v_vec = reverse_transform_matrix * transformed_v_vec;
    float scale_factor = v_vec.norm() / reverse_transformed_v_vec.norm();
    reverse_transformed_v_vec *= scale_factor;

    cout << "final result:\n" << reverse_transformed_v_vec << endl;
}