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
        for (unsigned i = 0; i < particle_number; ++i) {
            ensemble.push_back(Particle((i + 1) * 0.001, (i + 1) * 0.001, (i + 1) * 0.001, i + 1, i + 1, i + 1, 5, 5, 5, 1e-8));
        }

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->a_x = 0;
            particle->a_y = 0;
            particle->a_z = 0;
            for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
                // exclude the current particle
                if (&(*particle_ptr) != &(*particle)) {
                    particle->calculate_distance_value(*particle_ptr);
                    particle->acceleration(*particle_ptr);
                }
            }
            particle->former_a_x = particle->a_x;
            particle->former_a_y = particle->a_y;
            particle->former_a_z = particle->a_z;
        }
    }

    // use ensemble[index] to call particle's method
    Particle& operator[] (const int index) {
        return ensemble[index];
    }


    // basic interation step of velocity verlet
    void iteration(Particle& particle) {
        particle.movement();
        particle.a_x = 0;
        particle.a_y = 0;
        particle.a_z = 0;
        double distance_value = 0;
        particle.potential_value = 0;
        for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
            // exclude the current particle
            if (&(*particle_ptr) != &particle) {
                particle.calculate_distance_value(*particle_ptr);
                particle.acceleration(*particle_ptr);
            }
        }
        particle.velocity();
        particle.kinetic();

        particle.former_a_x = particle.a_x;
        particle.former_a_y = particle.a_y;
        particle.former_a_z = particle.a_z;

        ensemble_potential += particle.potential_value;
        ensemble_kinetic += particle.kinetic_value;
    }


    void output(ofstream& fout) {
        fout << ensemble_potential << "\t" << ensemble_kinetic << "\t" << (ensemble_potential + ensemble_kinetic) << "\n";
    }


    void execute(const unsigned time, const unsigned index, Box& box, ofstream& particle_out, ofstream& ensemble_out) {
        int i = 0;
        while (i <= time) {
            ensemble_kinetic = 0;
            ensemble_potential = 0;
            ensemble[index].output(particle_out);
            for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
                this->iteration(*particle_ptr);
                particle_ptr->rebounce(box);
            }
            output(ensemble_out);
            ++i;
        }
    }

private:
    
    vector<Particle> ensemble;

    unsigned particle_number;

    double ensemble_potential;

    double ensemble_kinetic;

};

#endif