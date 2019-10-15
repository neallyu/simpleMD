#ifndef PARTICLE_H
#define PARTICLE_H

#include <iostream>
#include <fstream>
#include <cmath>
#include <stdexcept>

using namespace std;

class Particle {

friend class Ensemble;
friend class Box;

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
            pow((pos_x - other.pos_x), 2) + 
            pow((pos_y - other.pos_y), 2) +
            pow((pos_z - other.pos_z), 2));
    }

    // calculate acceleration of the particle from interaction (particle version)
    void interact(const Particle& other) {
        calculate_distance_value(other);
        // a_x = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
        //     6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_x - pos_x) / mass;
        // a_y = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
        //     6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_y - pos_y) / mass;
        // a_z = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
        //     6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_z - pos_z) / mass;

        // another version of F(r) function to see if the speed is faster than before
        a_x = -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_x - pos_x) / mass;
        a_y = -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_y - pos_y) / mass;
        a_z = -4 * epsilon * (12 * pow(sigma / distance_value, 12) / distance_value - 
            6 * pow(sigma / distance_value, 6) / distance_value) * (other.pos_z - pos_z) / mass;
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

    // execute movement
    void movement() {
        // uniformly accelerated motion
        pos_x += v_x * time_interval + 0.5 * a_x * pow(time_interval, 2);
        pos_y += v_y * time_interval + 0.5 * a_y * pow(time_interval, 2);
        pos_z += v_z * time_interval + 0.5 * a_z * pow(time_interval, 2);
        v_x += a_x * time_interval;
        v_y += a_y * time_interval;
        v_z += a_z * time_interval;
    }

    // calculate the current kinetic energy of the particle
    void kinetic() {
        kinetic_value = 0.5 * mass * (pow(v_x, 2) + pow(v_y, 2) + pow(v_z, 2));
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

    double potential_value;

    double kinetic_value;

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

    double distance_value;

    double mass;

    double epsilon;

    double sigma;

    double time_interval;

};

#endif