#ifndef PARTICLE_H
#define PARTICLE_H

class Particle {

friend class Ensemble;
friend class CG;
friend double distance(Particle&, Particle&);
friend double distance2(Particle&, Particle&);
friend double calc_velocity(Particle&);

public:
    // Initializer, receive the initial status of the particle
    Particle(double _time_interval);

    // Copy initializer
    Particle(const Particle&);

    // calculate x from v and a
    void movement();

    // calculate velocity from a
    void velocity();

    // calculate the current kinetic energy of the particle
    inline void kinetic();

protected:
    // position
    double pos_x;
    double pos_y;
    double pos_z;

    // velocity
    double v_x;
    double v_y;
    double v_z;

    // acceleration of step A
    double a_x_A;
    double a_y_A;
    double a_z_A;

    // acceleration of step B 
    double a_x_B;
    double a_y_B;
    double a_z_B;

    double potential_value;
    double kinetic_value;

    const double time_interval;
};


#endif