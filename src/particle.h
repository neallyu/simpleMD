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
            mass(_mass), epsilon(_epsilon), sigma(_sigma), time_interval(_time_interval),
            sigma_6(_sigma * _sigma * _sigma * _sigma * _sigma * _sigma), 
            sigma_12(_sigma * _sigma * _sigma * _sigma * _sigma * _sigma * 
                _sigma * _sigma * _sigma * _sigma * _sigma * _sigma) {
                if (pos_x * pos_y * pos_z <= 0) {
                    throw runtime_error("Error: negative initial position");
                }
                if (mass <= 0) {
                    throw runtime_error("Error: invalid mass");
                }
                pos_x_A = pos_x - v_x * time_interval;
                pos_y_A = pos_y - v_y * time_interval;
                pos_z_A = pos_z - v_z * time_interval;
            }

    Particle(const Particle& other): v_x(other.v_x), v_y(other.v_y), v_z(other.v_z), pos_x(other.pos_x), pos_y(other.pos_y), pos_z(other.pos_z),
        mass(other.mass), epsilon(other.epsilon), sigma(other.sigma), time_interval(other.time_interval), sigma_6(other.sigma_6), sigma_12(other.sigma_12),
        pos_x_A(other.pos_x_A), pos_y_A(other.pos_y_A), pos_z_A(other.pos_z_A) { }

    // execute movement
    // void movement() {
    //     // Euler algorithm
    //     pos_x = pos_x + v_x * time_interval + 0.5 * a_x_A * time_interval * time_interval;
    //     pos_y = pos_y + v_y * time_interval + 0.5 * a_y_A * time_interval * time_interval;
    //     pos_z = pos_z + v_z * time_interval + 0.5 * a_z_A * time_interval * time_interval;
    // }


    void movement_new() {
        double x = 2 * pos_x - pos_x_A + a_x * time_interval * time_interval;
        double y = 2 * pos_y - pos_y_A + a_y * time_interval * time_interval;
        double z = 2 * pos_z - pos_z_A + a_z * time_interval * time_interval;

        v_x = (x - pos_x_A) / (2 * time_interval);
        v_y = (y - pos_y_A) / (2 * time_interval);
        v_z = (z - pos_z_A) / (2 * time_interval);

        kinetic();

        pos_x_A = pos_x;
        pos_y_A = pos_y;
        pos_z_A = pos_z;

        pos_x = x;
        pos_y = y;
        pos_z = z;
    }

    // void velocity() {                
    //     v_x = v_x + (a_x_A + a_x_B) * 0.5 * time_interval;
    //     v_y = v_y + (a_y_A + a_y_B) * 0.5 * time_interval;
    //     v_z = v_z + (a_z_A + a_z_B) * 0.5 * time_interval;
    // }

        // calculate the current kinetic energy of the particle
    void kinetic() {
        kinetic_value = 0.5 * mass * (v_x * v_x + v_y * v_y + v_z * v_z);
    }

    // calculate total energy (particle version)
    double total_energy() {
        return kinetic_value + potential_value;
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

    // output to file
    // void output(ofstream& fout) {
    //     fout << pos_x << "\t" << pos_y << "\t" << pos_z << "\t" << v_x << "\t" << v_y << "\t" << v_z << "\t" 
    //         << a_x_B << "\t" << a_y_B << "\t" << a_z_B << "\t" << potential_value << "\t" << kinetic_value << "\t" << total_energy() << "\n";
    // }
    void output(ofstream& fout) {
        fout << pos_x << "\t" << pos_y << "\t" << pos_z << "\t" << v_x << "\t" << v_y << "\t" << v_z << "\t" 
            << a_x << "\t" << a_y << "\t" << a_z << "\t" << potential_value << "\t" << kinetic_value << "\t" << total_energy() << "\n";
    }

    // print on terminal
    // void print() {
    //     cout << "position_x: " << pos_x << "\tv_x: " << v_x << "\ta_x: " << a_x_B << "\tposition_y: " << pos_y 
    //         << "\tv_y: " << v_y << "\ta_y: " << a_y_B << "\tposition_x: " << pos_z << "\tv_z: " << v_z << "\ta_z: " << a_z_B 
    //         << "\tpotential_value: " << potential_value << "\tkinetic value: " << kinetic_value << "\ttotal energy: " << total_energy() << endl;
    // }
    void print() {
        cout << "position_x: " << pos_x << "\tv_x: " << v_x << "\ta_x: " << a_x << "\tposition_y: " << pos_y 
            << "\tv_y: " << v_y << "\ta_y: " << a_y << "\tposition_x: " << pos_z << "\tv_z: " << v_z << "\ta_z: " << a_z 
            << "\tpotential_value: " << potential_value << "\tkinetic value: " << kinetic_value << "\ttotal energy: " << total_energy() << endl;
    }

protected:
    // position
    double pos_x;
    double pos_y;
    double pos_z;

    // position A
    double pos_x_A;
    double pos_y_A;
    double pos_z_A;

    // velocity
    double v_x;
    double v_y;
    double v_z;

    // acceleration
    double a_x;
    double a_y;
    double a_z;

    // // acceleration of step A
    // double a_x_A;
    // double a_y_A;
    // double a_z_A;

    // // acceleration of step B 
    // double a_x_B;
    // double a_y_B;
    // double a_z_B;

    double distance_value;
    double potential_value;
    double kinetic_value;

    const double sigma_6;
    const double sigma_12;

    double mass;
    double epsilon;
    double sigma;

    double time_interval;

};


#endif
