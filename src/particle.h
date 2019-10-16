#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>
#include "box.h"

using namespace std;

class Particle {

friend class Ensemble;

public:
    // Initializer, receive the initial status of the particle
    Particle(double _v_x, double _v_y, double _v_z, double _pos_x, double _pos_y, double _pos_z,
        double _mass, double _epsilon, double _sigma, double _time_interval): 
            v_x(_v_x), v_y(_v_y), v_z(_v_z), pos_x(_pos_x), pos_y(_pos_y), pos_z(_pos_z),
            mass(_mass), epsilon(_epsilon), sigma(_sigma), time_interval(_time_interval) {
                if (pos_x * pos_y * pos_z <= 0) {
                    throw runtime_error("Error: negative initial position");
                }
                if (mass <= 0) {
                    throw runtime_error("Error: invalid mass");
                }
            }

    Particle(const Particle& other): v_x(other.v_x), v_y(other.v_y), v_z(other.v_z), pos_x(other.pos_x), pos_y(other.pos_y), pos_z(other.pos_z),
        mass(other.mass), epsilon(other.epsilon), sigma(other.sigma), time_interval(other.time_interval) { }

    // calculate the currenct distance between two particle
    void calculate_distance_value(const Particle& other) {
        distance_value = sqrt(
            (pos_x - other.pos_x) * (pos_x - other.pos_x) + 
            (pos_y - other.pos_y) * (pos_y - other.pos_y) +
            (pos_z - other.pos_z) * (pos_z - other.pos_z)
        );
    }

    // execute movement
    void movement() {
        // uniformly accelerated motion
        pos_x += v_x * time_interval + 0.5 * a_x * time_interval * time_interval;
        pos_y += v_y * time_interval + 0.5 * a_y * time_interval * time_interval;
        pos_z += v_z * time_interval + 0.5 * a_z * time_interval * time_interval;
    }

    // calculate acceleration of the particle from interaction (particle version)
    // to save the resource, the potential of current partcile is calculated at the same time
    void acceleration(const Particle& other) {
        calculate_distance_value(other);
        a_x += -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_x - pos_x) / mass;
        a_y += -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_y - pos_y) / mass;
        a_z += -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_z - pos_z) / mass;

        potential_value += 4 * epsilon * ( pow(sigma / distance_value, 12) - pow(sigma / distance_value, 6) );
    }

    void velocity() {                
        v_x += (a_x + former_a_x) * 0.5 * time_interval;
        v_y += (a_y + former_a_y) * 0.5 * time_interval;
        v_z += (a_z + former_a_z) * 0.5 * time_interval;
    }

    // output to file
    void output(ofstream& fout) {
        fout << pos_x << "\t" << pos_y << "\t" << pos_z << "\t" << v_x << "\t" << v_y  << "\t" << v_z << "\t" 
            << a_x << "\t" << a_y << "\t" << a_z << "\t" << potential_value << "\t" << kinetic_value << "\t" << total_energy() <<"\n";
    }

    // print on terminal
    void print() {
        cout << "position_x: " << pos_x << "\tv_x: " << v_x << "\ta_x: " << a_x << "\tposition_y: " << pos_y 
            << "\tv_y: " << v_y << "\ta_y: " << a_y << "\tposition_x: " << pos_z << "\tv_z: " << v_z << "\ta_z: " << a_z 
            << "\tpotential_value: " << potential_value << "\tkinetic value: " << kinetic_value << "\ttotal energy: " << total_energy() << endl;
    }

    // rebounce if particle hits the wall of box (particle version)
    void rebounce(Box& box) {
        if ((pos_x >= box.dim1 && v_x > 0) || 
            (pos_x <= 0 && v_x < 0) ) {
            v_x = -v_x;
        }
        if ((pos_y >= box.dim2 && v_y > 0) || 
            (pos_y <= 0 && v_y < 0) ) {
            v_y = -v_y;
        }
        if ((pos_z >= box.dim3 && v_z > 0) || 
            (pos_z <= 0 && v_z < 0) ) {
            v_z = -v_z;
        }
    }

    // calculate the current kinetic energy of the particle
    void kinetic() {
        kinetic_value = 0.5 * mass * (v_x * v_x + v_y * v_y + v_z * v_z);
    }

    // calculate potential between two particles (particle version)
    // note: must call "interact" before calculate the potential
    void potential(const Particle& other) {
        potential_value = 4 * epsilon * ( pow(sigma / distance_value, 12) - pow(sigma / distance_value, 6) );
    }

    // calculate total energy (particle version)
    double total_energy() {
        return kinetic_value + potential_value;
    }

protected:
    // position
    double pos_x;
    double pos_y;
    double pos_z;

    // velocity
    double v_x;
    double v_y;
    double v_z;

    // acceleration
    double a_x;
    double a_y;
    double a_z;

    // former accerleration for velocity verlet
    double former_a_x;
    double former_a_y;
    double former_a_z;

    double distance_value;
    double potential_value;
    double kinetic_value;

    double mass;
    double epsilon;
    double sigma;

    double time_interval;

};


#endif