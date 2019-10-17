#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
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
        ensemble.push_back(Particle(1e-2, 1e-2, 1e-2, 1, 1, 1, 1, 1e-3, 1, 1e-4));
        ensemble.push_back(Particle(-1e-2, -1e-2, -1e-2, 2, 2, 2, 1, 1e-3, 1, 1e-4));

        // Initialize acceleartion in step A
        for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
            for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                calc_acceleration(*particle1, *particle2);
                particle1->a_x_A = particle1->a_x_B;
                particle1->a_y_A = particle1->a_y_B;
                particle1->a_z_A = particle1->a_z_B;
                particle2->a_x_A = particle2->a_x_B;
                particle2->a_y_A = particle2->a_y_B;
                particle2->a_z_A = particle2->a_z_B;
            }
        }
    }

    // use ensemble[index] to call particle's method
    Particle& operator[] (const int index) {
        return ensemble[index];
    }

    void calc_acceleration(Particle& particle1, Particle& particle2) {
        double distance_value = sqrt(
            (particle1.pos_x - particle2.pos_x) * (particle1.pos_x - particle2.pos_x) + 
            (particle1.pos_y - particle2.pos_y) * (particle1.pos_y - particle2.pos_y) +
            (particle1.pos_z - particle2.pos_z) * (particle1.pos_z - particle2.pos_z)
        );
        double distance_value_6 = distance_value * distance_value * distance_value * distance_value * distance_value * distance_value;
        double distance_value_12 = distance_value_6 * distance_value_6;
        double force_x = 48 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - particle1.sigma_6 / 0.5 * distance_value_6) * 
            (particle1.pos_x - particle2.pos_x) / (distance_value * distance_value);
        double force_y = 48 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - particle1.sigma_6 / 0.5 * distance_value_6) * 
            (particle1.pos_y - particle2.pos_y) / (distance_value * distance_value);
        double force_z = 48 * particle1.epsilon * (particle1.sigma_12 / distance_value_12 - particle1.sigma_6 / 0.5 * distance_value_6) * 
            (particle1.pos_z - particle2.pos_z) / (distance_value * distance_value);

        particle1.a_x_B += force_x / particle1.mass;
        particle1.a_y_B += force_y / particle1.mass;
        particle1.a_z_B += force_z / particle1.mass;

        particle2.a_x_B -= force_x / particle2.mass;
        particle2.a_y_B -= force_y / particle2.mass;
        particle2.a_z_B -= force_z / particle2.mass;

        particle1.potential_value += 4 * particle1.epsilon * ( particle1.sigma_12 / distance_value_12 - particle1.sigma_6 / distance_value_6 );
        particle2.potential_value += 4 * particle2.epsilon * ( particle2.sigma_12 / distance_value_12 - particle2.sigma_6 / distance_value_6 );
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
                particle->movement();

                // Initialize step B
                particle->a_x_B = 0;
                particle->a_y_B = 0;
                particle->a_z_B = 0;
                particle->potential_value = 0;
            }

            // calculate acceleration of step B
            for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
                for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                    calc_acceleration(*particle1, *particle2);
                }
            }

            // calculate velocity of step B
            for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
                particle->velocity();
                particle->kinetic();
                ensemble[index].output(particle_out);
                ensemble_potential += particle->potential_value;
                ensemble_kinetic += particle->kinetic_value;
                particle->rebounce(box);
            }
            output(ensemble_out);
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
