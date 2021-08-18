#include "particle.h"


Particle::Particle(double _time_interval): 
    v_x(0), v_y(0), v_z(0), pos_x(0), pos_y(0), pos_z(0), a_x_B(0), a_y_B(0), a_z_B(0), time_interval(_time_interval) { }

// Copy initializer
Particle::Particle(const Particle& other): v_x(other.v_x), v_y(other.v_y), v_z(other.v_z), pos_x(other.pos_x), pos_y(other.pos_y),
    pos_z(other.pos_z), a_x_B(other.a_x_B), a_y_B(other.a_y_B), a_z_B(other.a_z_B), time_interval(other.time_interval) { }


// Particle& Particle::operator=(const Particle& other) {
//     v_x = other.v_x;
//     v_y = other.v_y;
//     v_z = other.v_z;
//     pos_x = other.pos_x;
//     pos_y = other.pos_y;
//     pos_z = other.pos_z;
//     a_x_B = other.a_x_B;
//     a_y_B = other.a_y_B;
//     a_z_B = other.a_z_B;
//     time_interval = other.time_interval;
//     return *this;
// }


// execute movement
void Particle::movement() {
    // Use velocity verlet algorithm
    pos_x += v_x * time_interval + 0.5 * a_x_A * time_interval * time_interval;
    pos_y += v_y * time_interval + 0.5 * a_y_A * time_interval * time_interval;
    pos_z += v_z * time_interval + 0.5 * a_z_A * time_interval * time_interval;
}


void Particle::velocity() {
    // Use velocity verlet algorithm 
    v_x += (a_x_A + a_x_B) * 0.5 * time_interval;
    v_y += (a_y_A + a_y_B) * 0.5 * time_interval;
    v_z += (a_z_A + a_z_B) * 0.5 * time_interval;
}

// calculate the current kinetic energy of the particle
void Particle::kinetic() {
    kinetic_value = 0.5 * (v_x * v_x + v_y * v_y + v_z * v_z);
}
