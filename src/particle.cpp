#include "particle.h"

Particle::Particle(double _mass, double _epsilon, double _sigma, double _time_interval): 
    v_x(0), v_y(0), v_z(0), pos_x(0), pos_y(0), pos_z(0),
    mass(_mass), epsilon(_epsilon), sigma(_sigma), time_interval(_time_interval),
    rcut(2.5 * sigma), ecut(4.0 * (pow(sigma / rcut, 12) - pow(sigma / rcut, 6))),
    rlist2(sigma * sigma * 12.25) {
        if (mass <= 0) {
            throw runtime_error("Error: invalid mass");
        }
}

// Copy initializer
Particle::Particle(const Particle& other): v_x(other.v_x), v_y(other.v_y), v_z(other.v_z), pos_x(other.pos_x), 
    pos_y(other.pos_y), pos_z(other.pos_z),mass(other.mass), epsilon(other.epsilon), sigma(other.sigma),
    time_interval(other.time_interval), rcut(other.rcut), ecut(other.ecut), rlist2(other.rlist2) { }


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
    kinetic_value = 0.5 * mass * (v_x * v_x + v_y * v_y + v_z * v_z);
}

// output to file
void Particle::output(ofstream& fout) {
    fout << pos_x << "\t" << pos_y << "\t" << pos_z << "\t" << v_x << "\t" << v_y << "\t" << v_z << "\t" 
        << a_x_B << "\t" << a_y_B << "\t" << a_z_B << "\t" << potential_value << "\t" << kinetic_value 
        << "\t" << potential_value + kinetic_value << "\n";
}