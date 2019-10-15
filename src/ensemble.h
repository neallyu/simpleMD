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
    // Ensemble(const unsigned _particle_number): particle_number(_particle_number) {
    //     for (unsigned i = 0; i < particle_number; ++i) {
    //         ensemble.push_back(Particle((i + 1) * 0.001, (i + 1) * 0.001, (i + 1) * 0.001, i + 1, i + 1, i + 1, 5, 10, 2, 1e-5));
    //     }
    // }

    Ensemble(const unsigned _particle_number): particle_number(_particle_number) {
        for (unsigned i = 0; i < particle_number; ++i) {
            ensemble.push_back(Particle_Energy_Compensated((i + 1) * 0.001, (i + 1) * 0.001, (i + 1) * 0.001, i + 1, i + 1, i + 1, 5, 10, 2, 1e-5));
        }
    }

    // use ensemble[index] to call particle's method
    Particle& operator[] (const int index) {
        return ensemble[index];
    }

    // calculate the total acceleration from other particles (ensemble version)
    // to save the resource, the potential of current partcile is calculated at the same time
    void interact(Particle& particle) {
        particle.a_x = 0;
        particle.a_y = 0;
        particle.a_z = 0;
        double distance_value = 0;
        particle.potential_value = 0;
        for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
            // exclude the current particle
            if (&(*particle_ptr) != &particle) {
                particle.calculate_distance_value(*particle_ptr);
                particle.interact(*particle_ptr);

                // caluculate potential value
                particle.potential_value += 4 * particle.epsilon * ( pow(particle.sigma / particle.distance_value, 12) - pow(particle.sigma / particle.distance_value, 6) );
            }
        }
    }


    void interact(Particle_Energy_Compensated& particle) {
        particle.a_x = 0;
        particle.a_y = 0;
        particle.a_z = 0;
        double distance_value = 0;
        particle.potential_value = 0;
        for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
            // exclude the current particle
            if (&(*particle_ptr) != &particle) {
                particle.calculate_distance_value(*particle_ptr);
                particle.interact(*particle_ptr);

                // caluculate potential value
                particle.potential_value += 4 * particle.epsilon * ( pow(particle.sigma / particle.distance_value, 12) - pow(particle.sigma / particle.distance_value, 6) );
            }
        }
        particle.former_potential_value = particle.potential_value;
        particle.kinetic();
        particle.former_kinetic_value = particle.kinetic_value;
    }


    void energy_compensate(Particle_Energy_Compensated& particle) {
        particle.potential_value = 0;
        for (auto particle_ptr = ensemble.begin(); particle_ptr != ensemble.end(); ++particle_ptr) {
            // exclude the current particle
            if (&(*particle_ptr) != &particle) {
                particle.calculate_distance_value(*particle_ptr);
                particle.potential_value += 4 * particle.epsilon * ( pow(particle.sigma / particle.distance_value, 12) - pow(particle.sigma / particle.distance_value, 6) );
            }
        }
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
                interact(*particle_ptr);
                particle_ptr->movement();
                energy_compensate(*particle_ptr);
                *particle_ptr->energy_compensate();
                ensemble_potential += particle_ptr->potential_value;
                ensemble_kinetic += particle_ptr->kinetic_value;
                particle_ptr->rebounce(box);
            }
            output(ensemble_out);
            ++i;
        }
    }

private:
    
    vector<Particle_Energy_Compensated> ensemble;

    unsigned particle_number;

    double ensemble_potential;

    double ensemble_kinetic;

};

#endif