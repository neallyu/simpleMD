#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <cmath>
// #include <cstdlib>
// #include <ctime>
#include "box.h"
#include "particle.h"

using namespace std;

class Ensemble {

public:

    Ensemble(const unsigned _particle_number): particle_number(_particle_number) {
        // for (unsigned i = 0; i < particle_number; ++i) {
        //     ensemble.push_back(Particle((i + 1) * 1e-2, (i + 1) * 1e-2, (i + 1) * 1e-2, i + 1, i + 1, i + 1, 1, 1e-1, 1, 2e-6));
        // }
        ensemble.push_back(Particle(1e-3, 1e-3, 1e-3, 1, 1, 1, 1, 1e-3, 1, 1e-3));
        ensemble.push_back(Particle(-1e-3, -1e-3, -1e-3, 2, 2, 2, 1, 1e-3, 1, 1e-3));

        // Initialize acceleartion in step A
        // for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
        //     for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
        //         calc_acceleration(*particle1, *particle2);
        //         particle1->a_x_A = particle1->a_x_B;
        //         particle1->a_y_A = particle1->a_y_B;
        //         particle1->a_z_A = particle1->a_z_B;
        //         particle2->a_x_A = particle2->a_x_B;
        //         particle2->a_y_A = particle2->a_y_B;
        //         particle2->a_z_A = particle2->a_z_B;
        //     }
        // }
        for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
            for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                calc_acceleration_new(*particle1, *particle2);
            }
        }
    }

    // use ensemble[index] to call particle's method
    Particle& operator[] (const int index) {
        return ensemble[index];
    }

    // void calc_acceleration(Particle& particle1, Particle& particle2) {
    //     double distance_value = sqrt(
    //         (particle1.pos_x - particle2.pos_x) * (particle1.pos_x - particle2.pos_x) + 
    //         (particle1.pos_y - particle2.pos_y) * (particle1.pos_y - particle2.pos_y) +
    //         (particle1.pos_z - particle2.pos_z) * (particle1.pos_z - particle2.pos_z)
    //     );
    //     double distance_value_6 = distance_value * distance_value * distance_value * distance_value * distance_value * distance_value;
    //     double distance_value_12 = distance_value_6 * distance_value_6;
    //     double force_x = 48.0 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - 0.5 * particle1.sigma_6 / distance_value_6) * 
    //         (particle1.pos_x - particle2.pos_x) / (distance_value * distance_value);
    //     double force_y = 48.0 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - 0.5 * particle1.sigma_6 / distance_value_6) * 
    //         (particle1.pos_y - particle2.pos_y) / (distance_value * distance_value);
    //     double force_z = 48.0 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - 0.5 * particle1.sigma_6 / distance_value_6) * 
    //         (particle1.pos_z - particle2.pos_z) / (distance_value * distance_value);

    //     particle1.a_x_B += force_x / particle1.mass;
    //     particle1.a_y_B += force_y / particle1.mass;
    //     particle1.a_z_B += force_z / particle1.mass;

    //     particle2.a_x_B -= force_x / particle2.mass;
    //     particle2.a_y_B -= force_y / particle2.mass;
    //     particle2.a_z_B -= force_z / particle2.mass;

    //     particle1.potential_value += 4.0 * particle1.epsilon * ( particle1.sigma_12 / distance_value_12 - particle1.sigma_6 / distance_value_6 );
    //     particle2.potential_value += 4.0 * particle2.epsilon * ( particle2.sigma_12 / distance_value_12 - particle2.sigma_6 / distance_value_6 );
    // }

    void calc_acceleration_new(Particle &particle1, Particle &particle2) {

        double dx = particle1.pos_x - particle2.pos_x;
        double dy = particle1.pos_y - particle2.pos_y;
        double dz = particle1.pos_z - particle2.pos_z;

        dx -= 10 * round(dx / 10);
        dy -= 10 * round(dy / 10);
        dz -= 10 * round(dz / 10);

        double r2 = dx * dx + dy * dy + dz * dz;

        double r2i = particle1.sigma * particle1.sigma / r2;

        double r6i = pow(r2i, 3);

        double force_x = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dx / r2;
        double force_y = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dy / r2;
        double force_z = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dz / r2;


        particle1.a_x += force_x / particle1.mass;
        particle1.a_y += force_y / particle1.mass;
        particle1.a_z += force_z / particle1.mass;

        particle2.a_x -= force_x / particle2.mass;
        particle2.a_y -= force_y / particle2.mass;
        particle2.a_z -= force_z / particle2.mass;

        particle1.potential_value += 4.0 * particle1.epsilon * r6i * ( r6i - 1 );
        particle2.potential_value += 4.0 * particle2.epsilon * r6i * ( r6i - 1 );
    }


    // basic interation step of velocity verlet
    void iteration(const unsigned time, const unsigned index, Box& box, ofstream& particle_out, ofstream& ensemble_out) {
        int i = 0;
        while (i <= time) {

            // Initialize ensemble energy
            ensemble_kinetic = 0;
            ensemble_potential = 0;

            // Execute movement of step A & Initialize acceleration of step B
            for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
                // step A
                particle->movement_new();

                // Initialize step B
                particle->a_x = 0;
                particle->a_y = 0;
                particle->a_z = 0;
                particle->potential_value = 0;
            }

            // calculate acceleration of step B
            for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
                for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                    calc_acceleration_new(*particle1, *particle2);
                }
            }

            // calculate velocity of step B
            for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
                // particle->velocity();
                // particle->kinetic();
                // if (i % 1000 == 0) {
                ensemble[index].output(particle_out);
                // }
                ensemble_potential += particle->potential_value;
                ensemble_kinetic += particle->kinetic_value;
            }
            // if (i % 1000 == 0) {
            output(ensemble_out);
            // }
            ++i;
        }
    }


    void output(ofstream& fout) {
        fout << ensemble_potential << "\t" << ensemble_kinetic << "\t" << (ensemble_potential + ensemble_kinetic) << "\n";
    }

private:
    
    vector<Particle> ensemble;

    unsigned particle_number;

    double ensemble_potential;

    double ensemble_kinetic;

};

#endif
