#ifndef PARTICLE_H
#define PARTICLE_H

#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace std;

class Particle {

friend class Ensemble;
friend class Neighborlist;

public:
    // Initializer, receive the initial status of the particle
    Particle(double, double, double, double);

    // Copy initializer
    Particle(const Particle&);

    // initialize lattice postion of particle
    inline void lattice_pos(int); 
    
    // calculate x from v and a
    inline void movement();

    // calculate velocity from a
    inline void velocity();

    // calculate the current kinetic energy of the particle
    inline void kinetic();

    // output to file
    inline void output(ofstream&);

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

    const double mass;
    const double epsilon;
    const double sigma;

    // cutoff distance and corresponding potential energy
    const double rcut;
    const double ecut;

    // distance threshold of neighborlist
    const double rlist2;

    const double time_interval;
};

#endif
