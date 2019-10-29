#ifndef PARTICLE_H
#define PARTICLE_H

#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace std;

class Particle {

friend class Ensemble;
friend class Neighborlist;
friend double distance(Particle&, Particle&, double);
friend double distance2(Particle&, Particle&, double);

public:
    // Initializer, receive the initial status of the particle
    Particle(double _time_interval);

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


Particle::Particle(double _time_interval): 
    v_x(0), v_y(0), v_z(0), pos_x(0), pos_y(0), pos_z(0), time_interval(_time_interval) { }

// Copy initializer
Particle::Particle(const Particle& other): v_x(other.v_x), v_y(other.v_y), v_z(other.v_z), pos_x(other.pos_x), pos_y(other.pos_y),
    pos_z(other.pos_z), time_interval(other.time_interval) { }


void Particle::lattice_pos(int n) {
    int i;
    for (i = 0; i * i * i < n; ++i);
    int j = n - (i - 1)*(i - 1)*(i - 1);
    if (j <= i * i) {
        pos_x = i - 1;
        pos_y = (j - 1) / i;
        pos_z = (j - 1) % i;
    }
    else if (j <= i * (2 * i - 1)) {
        pos_x = (j - i * i - 1) / i;
        pos_y = i - 1;
        pos_z = (j - i * i - 1) % i;
    }
    else {
        pos_x = (j - i * (2 * i - 1) - 1) / (i - 1);
        pos_y = (j - i * (2 * i - 1) - 1) % (i - 1);
        pos_z = i - 1;
    }
}
    
    
// execute movement
void Particle::movement() {
    // Euler algorithm
    pos_x = pos_x + v_x * time_interval + 0.5 * a_x_A * time_interval * time_interval;
    pos_y = pos_y + v_y * time_interval + 0.5 * a_y_A * time_interval * time_interval;
    pos_z = pos_z + v_z * time_interval + 0.5 * a_z_A * time_interval * time_interval;
}


void Particle::velocity() {                
    v_x = v_x + (a_x_A + a_x_B) * 0.5 * time_interval;
    v_y = v_y + (a_y_A + a_y_B) * 0.5 * time_interval;
    v_z = v_z + (a_z_A + a_z_B) * 0.5 * time_interval;
}

// calculate the current kinetic energy of the particle
void Particle::kinetic() {
    kinetic_value = 0.5 * (v_x * v_x + v_y * v_y + v_z * v_z);
}


#endif
