#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <cmath>
#include <random>
#include "box.h"
#include "particle.h"
#include "scale.h"

using namespace std;

class Ensemble {

public:

    Ensemble(const unsigned _particle_number, double temp, double time_interval, Box& box): particle_number(_particle_number), TEMP(temp) {
        for (unsigned i = 0; i < particle_number; ++i) {
            ensemble.push_back(Particle(0, 0, 0, 0, 0, 0, amu(32), 407.6, nm(0.302), time_interval));
        }

        init_particle(box);

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
        for (auto particle1 = ensemble.begin(); particle1 != ensemble.end(); ++particle1) {
            for (auto particle2 = particle1 + 1; particle2 != ensemble.end(); ++particle2) {
                calc_acceleration(*particle1, *particle2);
            }
        }
    }

    // use ensemble[index] to call particle's method
    Particle& operator[] (const int index) {
        return ensemble[index];
    }

    void calc_acceleration(Particle& particle1, Particle& particle2) {

        double dx = particle1.pos_x - particle2.pos_x;
        double dy = particle1.pos_y - particle2.pos_y;
        double dz = particle1.pos_z - particle2.pos_z;

        dx -= 10 * round(dx / 10);
        dy -= 10 * round(dy / 10);
        dz -= 10 * round(dz / 10);

        double r2 = dx * dx + dy * dy + dz * dz;

        if (r2 <= particle1.rcut * particle1.rcut) {

            double r2i = particle1.sigma * particle1.sigma / r2;

            double r6i = pow(r2i, 3);

            double force_x = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dx / r2;
            double force_y = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dy / r2;
            double force_z = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dz / r2;


            particle1.a_x_B += force_x / particle1.mass;
            particle1.a_y_B += force_y / particle1.mass;
            particle1.a_z_B += force_z / particle1.mass;

            particle2.a_x_B -= force_x / particle2.mass;
            particle2.a_y_B -= force_y / particle2.mass;
            particle2.a_z_B -= force_z / particle2.mass;

            particle1.potential_value += 4.0 * particle1.epsilon * r6i * ( r6i - 1 ) - particle1.ecut;
            particle2.potential_value += 4.0 * particle2.epsilon * r6i * ( r6i - 1 ) - particle2.ecut;
        }

    }

    void init_particle(Box &box) {
        default_random_engine random_generator;
        uniform_real_distribution<double> displacement(0.0, 1.0);  //distribution generator

        double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0);
        double sumv2(0.0);

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->pos_x = displacement(random_generator) * box.dim1;
            cout << "Initial position: " << particle->pos_x << "\t";
            particle->pos_y = displacement(random_generator) * box.dim2;
            cout << particle->pos_y << "\t";
            particle->pos_z = displacement(random_generator) * box.dim3;
            cout << particle->pos_z << endl;

            particle->v_x = nm(displacement(random_generator) - 0.5);
            particle->v_y = nm(displacement(random_generator) - 0.5);
            particle->v_z = nm(displacement(random_generator) - 0.5);

            sumv_x += particle->v_x;
            sumv_y += particle->v_y;
            sumv_z += particle->v_z;

            sumv2 += sumv_x * sumv_x + sumv_y * sumv_y + sumv_z * sumv_z;
        }

        sumv_x /= particle_number;
        sumv_y /= particle_number;
        sumv_z /= particle_number;

        sumv2 /= particle_number;
        
        double fs = sqrt(3 * kb * TEMP / (sumv2 * ensemble[0].mass));

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->v_x = (particle->v_x - sumv_x) * fs;
            cout << "Initial velocity: " << particle->v_x << "\t";
            particle->v_y = (particle->v_y - sumv_y) * fs;
            cout << particle->v_y << "\t";
            particle->v_z = (particle->v_z - sumv_z) * fs;
            cout << particle->v_z << endl;
        }
    }

    // void calc_acceleration_new(Particle &particle1, Particle &particle2) {

    //     double dx = particle1.pos_x - particle2.pos_x;
    //     double dy = particle1.pos_y - particle2.pos_y;
    //     double dz = particle1.pos_z - particle2.pos_z;

    //     dx -= 10 * round(dx / 10);
    //     dy -= 10 * round(dy / 10);
    //     dz -= 10 * round(dz / 10);

    //     double r2 = dx * dx + dy * dy + dz * dz;

    //     double r2i = particle1.sigma * particle1.sigma / r2;

    //     double r6i = pow(r2i, 3);

    //     double force_x = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dx / r2;
    //     double force_y = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dy / r2;
    //     double force_z = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dz / r2;


    //     particle1.a_x += force_x / particle1.mass;
    //     particle1.a_y += force_y / particle1.mass;
    //     particle1.a_z += force_z / particle1.mass;

    //     particle2.a_x -= force_x / particle2.mass;
    //     particle2.a_y -= force_y / particle2.mass;
    //     particle2.a_z -= force_z / particle2.mass;

    //     particle1.potential_value += 4.0 * particle1.epsilon * r6i * ( r6i - 1 );
    //     particle2.potential_value += 4.0 * particle2.epsilon * r6i * ( r6i - 1 );
    // }


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
                // if (i % 1000 == 0) {
                ensemble[index].output(particle_out);
                // }
                ensemble_potential += particle->potential_value;
                ensemble_kinetic += particle->kinetic_value;
                particle->rebounce(box);
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

    const double TEMP;


};

#endif
