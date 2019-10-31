#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include "particle.hpp"
#include "neighborlist.hpp"
#include "unit_conversion.hpp"
#include "radial_distribution_function.hpp"
#include "utils.hpp"

using namespace std;

class Ensemble {

friend class Neighborlist;
friend class Rdf;

public:
    // box(angstrom), temp(K), sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
    Ensemble(const unsigned _particle_number, double sigma, double epsilon, double mass, double temp, double time_interval, 
    unsigned long _TIME, double box);

    ~Ensemble();

    // return particle of that index
    inline Particle& operator[] (const int);

    // return distance^2 between two particles
    inline double calc_distance2(Particle&, Particle&);

    // calculation <a> between a pair of particles
    inline void calc_acceleration(Particle&, Particle&);
    inline void iteration();
    inline void energy_output(unsigned long i, ofstream& fout);
    inline void particle_movement_output(unsigned long i, Particle& particle, ofstream& fout);


private:
    // Initialize unit conversion
    Unit_conversion unit;

    // reduced unit (calculated from sigma, epsilon, mass)
    const double TEMP;
    const double BOX;
    const double TIME_INTERVAL;

    // main container of the particle ensemble
    unsigned particle_number;
    vector<Particle> ensemble;

    // cutoff distance, cutoff potential energy and distance threshold of neighborlist
    const double rcut;
    const double ecut;
    const double rlist2;

    Neighborlist nlist;
    bool need_update_nlist;

    Rdf rdf;

    const double TIME;
    const unsigned long SAMPLE_RATE;

    double ensemble_potential;
    double ensemble_kinetic;
    ofstream ensemble_out;
    ofstream particle_out;
};

// box(angstrom), temp(K), sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
Ensemble::Ensemble(const unsigned _particle_number, double sigma, double epsilon, double mass, double temp, 
    double time_interval, unsigned long _TIME, double box): 
    particle_number(_particle_number), unit(sigma, epsilon, mass), TEMP(unit.reduced_temperature(temp)), 
    BOX(unit.reduced_distance(box * 1e-10)), TIME_INTERVAL(unit.reduced_time(time_interval * 1e-12)), 
    ensemble(_particle_number, Particle(TIME_INTERVAL)),rcut(2.5), ecut(4.0 * (pow(1 / rcut, 12) - pow(1 / rcut, 6))),
    rlist2(12.25), nlist(ensemble, BOX, rlist2), rdf(1000, BOX), TIME(unit.reduced_time(_TIME * 1e-9)), SAMPLE_RATE(TIME / 1000), 
    ensemble_out("../output/energy.csv"), particle_out("../output/particle.csv") {

        cout << "TIME_Interval: " << TIME_INTERVAL << endl;
        cout << "time: " << TIME << endl;
        cout << "sample rate: " << SAMPLE_RATE << endl;

        cout << "[MD LOG] " << get_current_time() << "\tEnsemble energy data output to \"../output/energy.csv\" ..." << endl;
        cout << "[MD LOG] " << get_current_time() << "\tParticle trajectory data output to \"../output/particle.csv\" ..." << endl;

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

        double fs = sqrt(3 * TEMP / sumv2);

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

Ensemble::~Ensemble() {
    ensemble_out.close();
    particle_out.close();

    cout << "[MD LOG] " << get_current_time() << "\tOutput file saved" << endl;
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

    if (r2 <= rcut * rcut) {
        double r2i = 1 / r2;
        double r6i = pow(r2i, 3);

        double force_x = 48.0 * r6i * r2i * (r6i - 0.5) * dx;
        double force_y = 48.0 * r6i * r2i * (r6i - 0.5) * dy;
        double force_z = 48.0 * r6i * r2i * (r6i - 0.5) * dz;

        particle1.a_x_B += force_x;
        particle1.a_y_B += force_y;
        particle1.a_z_B += force_z;

        particle2.a_x_B -= force_x;
        particle2.a_y_B -= force_y;
        particle2.a_z_B -= force_z;

        particle1.potential_value += (4.0 * r6i * ( r6i - 1 ) - ecut);
    }

    if (r2 > nlist.rlist2) {
        need_update_nlist = true;
    }
}


void Ensemble::energy_output(unsigned long i, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << unit.real_energy(ensemble_potential) << "    " << unit.real_energy(ensemble_kinetic) << "    " 
    << unit.real_energy(ensemble_potential + ensemble_kinetic) << endl;
}


// // output to file
void Ensemble::particle_movement_output(unsigned long i, Particle& particle, ofstream& fout) {
    fout << i * TIME_INTERVAL << "    " << unit.real_distance(particle.pos_x) << "    " << unit.real_distance(particle.pos_y) << "    " << unit.real_distance(particle.pos_z) << "    " 
        << unit.real_velocity(particle.v_x) << "    " << unit.real_velocity(particle.v_y) << "    " << unit.real_velocity(particle.v_z) << "    " 
        << unit.real_acceleration(particle.a_x_B) << "    " << unit.real_acceleration(particle.a_y_B) << "    " << unit.real_acceleration(particle.a_z_B) << "    " 
        << unit.real_energy(particle.potential_value) << "    " << unit.real_energy(particle.kinetic_value) 
        << "    " << unit.real_energy(particle.potential_value + particle.kinetic_value) << endl;
}


// basic interation step of velocity verlet
void Ensemble::iteration() {
    unsigned long i = 0;
    while (i <= TIME) {
        // Initialize ensemble energy
        ensemble_kinetic = 0;
        ensemble_potential = 0;

        // calculate acceleration of step B in neighbor list
        for (int i = 0; i < nlist.nlist.size(); ++i) {
            for (auto j = nlist.nlist[i].begin(); j != nlist.nlist[i].end(); ++j) {
                calc_acceleration(ensemble[i], ensemble[*j]);
            }
        }

        if (need_update_nlist == true) {
            nlist.update_neighbor_list(ensemble);
            need_update_nlist = false;
        }

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
        // output
        if (i % SAMPLE_RATE == 0) {
            // show the trajectory of one particle
            particle_movement_output(i, ensemble[1], particle_out);
            energy_output(i, ensemble_out);
            rdf.sample(ensemble);
        }
        ++i;
    }

    rdf.normalize(particle_number);
    rdf.output();
}


#endif
