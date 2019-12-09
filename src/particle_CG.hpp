#ifndef PARTICLE_CG_H
#define PARTICLE_CG_H

class Particle_CG {

friend class Ensemble;
friend double distance(Particle_CG&, Particle_CG&);

public:
    // Initializing all variant
    Particle_CG();

    // calculate the gradient with other particle
    double calc_gradient(Particle_CG&);

protected:
    // position of step A
    double pos_x_A;
    double pos_y_A;
    double pos_z_A;

    // position of step B
    double pos_x_B;
    double pos_y_B;
    double pos_z_B;

    // position of step C
    double pos_x_C;
    double pos_y_C;
    double pos_z_C;

    // direction of step A
    double d_x_A;
    double d_y_A;
    double d_z_A;

    // direction of step B
    double d_x_B;
    double d_y_B;
    double d_z_B;

    // direction of step C
    double d_x_C;
    double d_y_C;
    double d_z_C;

    // gradient of step A
    double g_x_A;
    double g_y_A;
    double g_z_A;

    // gradient of step B
    double g_x_B;
    double g_y_B;
    double g_z_B;

    // gradient of step C
    double g_x_C;
    double g_y_C;
    double g_z_C;
};


Particle_CG::Particle_CG(): 
    pos_x_A(0), pos_y_A(0), pos_z_A(0), pos_x_B(0), pos_y_B(0), pos_z_B(0), pos_x_C(0), pos_y_C(0), pos_z_C(0),
    d_x_A(0), d_y_A(0), d_z_A(0), d_x_B(0), d_y_B(0), d_z_B(0), d_x_C(0), d_y_C(0), d_z_C(0),
    g_x_A(0), g_y_A(0), g_z_A(0), g_x_B(0), g_y_B(0), g_z_B(0), g_x_C(0), g_y_C(0), g_z_C(0)
    { }


#endif
