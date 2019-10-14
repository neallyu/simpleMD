#ifndef PARTICLE_H
#define PARTICLE_H

#include <fstream>
#include <vector>
#include <cmath>
#include <stdexcept>
#include "integral.h"

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

    // check if two particles are the same
    // bool operator==(const Particle &other) {
    //     return pos_x == other.pos_x && pos_y == other.pos_y && pos_z == other.pos_z &&
    //             v_x == other.v_x && v_y == other.v_y && v_z == other.v_z &&
    //             a_x == other.a_x && a_y == other.a_y && a_z == other.a_z;
    // }

    // calculate acceleration of the particle from interaction
    void interact(const Particle& other) {
        double distance_value = sqrt(
            pow((pos_x - other.pos_x), 2) + 
            pow((pos_y - other.pos_y), 2) +
            pow((pos_z - other.pos_z), 2));
        a_x = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
            6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_x - pos_x) / mass;
        a_y = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
            6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_y - pos_y) / mass;
        a_z = -4 * epsilon * (12 * pow(sigma, 12) / pow(distance_value, 13) - 
            6 * pow(sigma, 6) / pow(distance_value, 7)) * (other.pos_z - pos_z) / mass;
    }

    // output to file
    void output(ofstream& fout) {
        fout << pos_x << "\t" << pos_y << "\t" << pos_z << "\t" << v_x << "\t" << v_y 
            << "\t" << v_z << "\t" << a_x << "\t" << a_y << "\t" << a_z <<"\n";
    }

    // print on terminal
    void print() {
        cout << "position_x: " << pos_x << "\tv_x: " << v_x << "\ta_x: " << a_x << "\tposition_y: " << pos_y << 
            "\tv_y: " << v_y << "\ta_y: " << a_y << "\tposition_x: " << pos_z << "\tv_z: " << v_z << "\ta_z: " << a_z << endl;
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

protected:

    // velocity
    double v_x;
    double v_y;
    double v_z;

    // position
    double pos_x;
    double pos_y;
    double pos_z;


    // acceleration
    double a_x;
    double a_y;
    double a_z;

    double mass;

    double epsilon;

    double sigma;

    double time_interval;

};



// class Particle_Energy_Corrected: public Particle {
// public:

//     // calculate potential between two particles
//     double potential(const double &_pos_x, const double &_pos_y, const double &_pos_z, 
//         const double &_other_pos_x, const double &_other_pos_y, const double &_other_pos_z) {

//         double distance_value = sqrt(
//             pow((_pos_x - _other_pos_x), 2) + 
//             pow((_pos_y - _other_pos_y), 2) +
//             pow((_pos_z - _other_pos_z), 2));
        
//         if (distance_value <= sigma) {
//             // calculate the intergration in (r, sigma)
//             return integral(epsilon, sigma, distance_value, sigma);
//         } else {
//             return -1 * integral(epsilon, sigma, sigma, distance_value);
//         }
//     }

//     // calculate kinetic energy of the particle
//     double kinetic(const double &_v_x, const double &_v_y, const double &_v_z) {
//         return 0.5 * mass * sqrt( pow(_v_x, 2) + pow(_v_y, 2) + pow(_v_z, 2) );
//     }

//     // calculate the velocity according to conservation of energy 
//     void interact(const Particle_Energy_Corrected& other) {
//         double velocity = sqrt(
//             (kinetic(former_v_x, former_v_y, former_v_z) + 
//             potential(former_pos_x, former_pos_y, former_pos_z, other.former_pos_x, other.former_pos_y, other.former_pos_z) -
//             potential(pos_x, pos_y, pos_z, other.pos_x, other.pos_y, other.pos_z) ) * 2 / mass
//         );

//         double velocity_to_be_corrected = sqrt(pow(v_x, 2) + pow(v_y, 2) + pow(v_z, 2));
//         v_x = velocity * v_x / velocity_to_be_corrected;
//         v_y = velocity * v_y / velocity_to_be_corrected;
//         v_z = velocity * v_z / velocity_to_be_corrected;
//     }

//     // calculate movement in given time granularity
//     void movement() {
//         // record former position and velocity
//         former_pos_x = pos_x;
//         former_pos_y = pos_y;
//         former_pos_z = pos_z;

//         former_v_x = v_x;
//         former_v_y = v_y;
//         former_v_z = v_z;

//         // uniformly accelerated motion
//         pos_x += v_x * time_interval + 0.5 * a_x * pow(time_interval, 2);
//         pos_y += v_y * time_interval + 0.5 * a_y * pow(time_interval, 2);
//         pos_z += v_z * time_interval + 0.5 * a_z * pow(time_interval, 2);
//         v_x += a_x * time_interval;
//         v_y += a_y * time_interval;
//         v_z += a_z * time_interval;
//     }

// private:
//     // former position
//     double former_pos_x;
//     double former_pos_y;
//     double former_pos_z;

//     //former velocity
//     double former_v_x;
//     double former_v_y;
//     double former_v_z;
// };

#endif