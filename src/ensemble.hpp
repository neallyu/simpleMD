#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include "particle.hpp"
#include "neighborlist.hpp"
#include "radial_distribution_function.hpp"

using namespace std;

class Ensemble {

friend class Neighborlist;
friend class Rdf;

public:
    // Initializer
    Ensemble(const unsigned, double, double, double);

    // return particle of that index
    inline Particle& operator[] (const int);

    // return distance^2 between two particles
    inline double calc_distance2(Particle&, Particle&);

    // calculation <a> between a pair of particles
    inline void calc_acceleration(Particle&, Particle&);

    inline void iteration(const unsigned, const unsigned, ofstream&, ofstream&);

    inline void output(ofstream&);


private:

    vector<Particle> ensemble;
    Neighborlist nlist;
    bool need_update_nlist;
    unsigned particle_number;
    double ensemble_potential;
    double ensemble_kinetic;
    const double TEMP;
    const double BOX;
    Rdf rdf;
};


Ensemble::Ensemble(const unsigned _particle_number, double temp, double time_interval, double box): 
    particle_number(_particle_number), TEMP(temp), BOX(box), 
    ensemble(_particle_number, Particle(32, 1, 1, time_interval)), 
    nlist(ensemble, box, ensemble[0].rlist2), rdf(1000, box) {

        default_random_engine random_generator;
        uniform_real_distribution<double> displacement(0.0, 1.0);  //distribution generator

        double sumv_x(0.0), sumv_y(0.0), sumv_z(0.0);
        double sumv2(0.0);

        int i = 0; // for lattice pos

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {

            particle->lattice_pos(i + 1);
            ++i;

            particle->pos_x += 0.01 * (displacement(random_generator) - 0.5);
            // cout << "Initial position: " << particle->pos_x << "\t";
            particle->pos_y += 0.01 * (displacement(random_generator) - 0.5);
            // cout << particle->pos_y << "\t";
            particle->pos_z += 0.01 * (displacement(random_generator) - 0.5);
            // cout << particle->pos_z << endl;

            particle->v_x = displacement(random_generator) - 0.5;
            particle->v_y = displacement(random_generator) - 0.5;
            particle->v_z = displacement(random_generator) - 0.5;

            sumv_x += particle->v_x;
            sumv_y += particle->v_y;
            sumv_z += particle->v_z;

            sumv2 += sumv_x * sumv_x + sumv_y * sumv_y + sumv_z * sumv_z;
        }

        sumv_x /= particle_number;
        sumv_y /= particle_number;
        sumv_z /= particle_number;

        sumv2 /= particle_number;

        double fs = sqrt(3 * TEMP / (sumv2 * ensemble[0].mass));

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            particle->v_x = (particle->v_x - sumv_x) * fs;
            // cout << "Initial velocity: " << particle->v_x << "\t";
            particle->v_y = (particle->v_y - sumv_y) * fs;
            // cout << particle->v_y << "\t";
            particle->v_z = (particle->v_z - sumv_z) * fs;
            // cout << particle->v_z << endl;
        }

        // Initialize neighborlist
        nlist.update_neighbor_list(ensemble);
        need_update_nlist = false;

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

        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            // execute x and v propagation
            particle->movement();
            particle->velocity();

            // record a_A and initialize a_B
            particle->a_x_A = particle->a_x_B;
            particle->a_y_A = particle->a_y_B;
            particle->a_z_A = particle->a_z_B;
            particle->a_x_B = 0;
            particle->a_y_B = 0;
            particle->a_z_B = 0;
            particle->potential_value = 0;
        }
}

// use ensemble[index] to call particle's method
Particle& Ensemble::operator[] (const int index) {
    return ensemble[index];
}


double Ensemble::calc_distance2(Particle& particle1, Particle& particle2) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    return dx * dx + dy * dy + dz * dz;
}


void Ensemble::calc_acceleration(Particle& particle1, Particle& particle2) {
    double dx = particle1.pos_x - particle2.pos_x;
    double dy = particle1.pos_y - particle2.pos_y;
    double dz = particle1.pos_z - particle2.pos_z;

    // periodic boundary conditions
    dx -= BOX * round(dx / BOX);
    dy -= BOX * round(dy / BOX);
    dz -= BOX * round(dz / BOX);

    double r2 = dx * dx + dy * dy + dz * dz;

    if (r2 <= particle1.rcut * particle1.rcut) {

        double r2i = particle1.sigma * particle1.sigma / r2;

        double r6i = pow(r2i, 3);

        double force_x = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dx / r2;
        double force_y = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dy / r2;
        double force_z = 48.0 * particle1.epsilon * r6i * (r6i - 0.5) * dz / r2;


        particle1.a_x_B += (force_x / particle1.mass);
        particle1.a_y_B += (force_y / particle1.mass);
        particle1.a_z_B += (force_z / particle1.mass);

        particle2.a_x_B -= (force_x / particle2.mass);
        particle2.a_y_B -= (force_y / particle2.mass);
        particle2.a_z_B -= (force_z / particle2.mass);

        particle1.potential_value += (4.0 * particle1.epsilon * r6i * ( r6i - 1 ) - particle1.ecut);
        particle2.potential_value += (4.0 * particle2.epsilon * r6i * ( r6i - 1 ) - particle2.ecut);
    }

    if (r2 > nlist.rlist2) {
        need_update_nlist = true;
    }
}


void Ensemble::output(ofstream& fout) {
    fout << ensemble_potential << "\t" << ensemble_kinetic << "\t" << (ensemble_potential + ensemble_kinetic) << "\n";
}


// basic interation step of velocity verlet
void Ensemble::iteration(const unsigned time, const unsigned index, ofstream& particle_out, ofstream& ensemble_out) {
    int i = 0;
    while (i <= time) {
        // Initialize ensemble energy
        ensemble_kinetic = 0;
        ensemble_potential = 0;

        // calculate acceleration of step B in neighbor list
        // #pragma omp parallel
        for (int i = 0; i < nlist.nlist.size(); ++i) {
            for (auto j = nlist.nlist[i].begin(); j != nlist.nlist[i].end(); ++j) {
                calc_acceleration(ensemble[i], ensemble[*j]);
            }
        }

        if (need_update_nlist == true) {
            nlist.update_neighbor_list(ensemble);
            need_update_nlist = false;
        }

        // #pragma omp parallel
        // calculate velocity of step B
        for (auto particle = ensemble.begin(); particle != ensemble.end(); ++particle) {
            // record energy
            particle->kinetic();
            ensemble_potential += particle->potential_value;
            ensemble_kinetic += particle->kinetic_value;

            // execute x and v propagation
            particle->movement();
            particle->velocity();

            // record a_A and initialize a_B
            particle->a_x_A = particle->a_x_B;
            particle->a_y_A = particle->a_y_B;
            particle->a_z_A = particle->a_z_B;
            particle->a_x_B = 0;
            particle->a_y_B = 0;
            particle->a_z_B = 0;
            particle->potential_value = 0;
        }

        // This is necessary because the potential energy of individual particle is calculated twice
        ensemble_potential /= 2;

        // output
        if (i % 100 == 0) {
            ensemble[index].output(particle_out);
            output(ensemble_out);

            rdf.sample(ensemble);
        }
        ++i;
    }

    rdf.normalize(particle_number);
    rdf.output();
}


#endif
