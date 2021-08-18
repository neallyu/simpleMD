#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include "particle.h"
#include "neighborlist.h"
#include "radial_distribution_function.h"
#include "property.h"
#include "utils.h"

using namespace std;

class Ensemble {

friend class Neighborlist;
friend class Rdf;
friend class Property;

public:
    // reduced unit
    Ensemble(const unsigned _particle_number, double init_temp, double set_temp, double time_interval, 
    double equilibration_time, double total_time, double rho, string _output_path);

    // close file stream
    ~Ensemble();

    // return object of particle of that index
    inline Particle& operator[] (const int index);

    // lattice position
    inline void lattice_pos();

    // calculation acceleration between a pair of particles
    inline void calc_acceleration(Particle& particle1, Particle& particle2);

    // correct the instaneous temperature by rescaling
    inline void rescale_temperature(double target_TEMP);

    inline void Andersen_thermostat(double collision_frequency);

    // main iteration
    void iteration();

    inline void energy_output(unsigned long i, ofstream& fout);
    inline void particle_movement_output(unsigned long i, Particle& particle, ofstream& fout);
    inline void temperature_output(unsigned long i, ofstream& fout);
    inline void msd_output(unsigned long i, double _MSD, ofstream& fout);

private:
    // all variable below are in reduced unit
    const double INIT_TEMP;                         // initial temperature before equilibration
    const double SET_TEMP;                          // set temperature after equilibration
    double TEMP;                                    // instaneous temperature
    unsigned particle_number;                       // particle number
    const double rho;                               // density
    const double BOX;                               // box dimension size
    const double TIME_INTERVAL;                     // machine time interval
    const double EQUILIBRATION_TIME;                // euquilibration time
    const double TOTAL_TIME;                        // total time
    const unsigned long EQUILIBRATION_ITERATION;    // iteration cycles of equilibration
    const unsigned long ITERATION;                  // total iteration cycles
    const unsigned long SAMPLE_RATE;                // sample rate defined as (iteration / 1000) such that the result contains 1000 points
    float ITERATION_PERCENTAGE;                      // percentage of main iteration
    vector<Particle> ensemble;                      // main container of the particle ensemble
    const double rcut;                              // cutoff distance defined as 2.5 (reduced unit)
    const double ecut;                              // cutoff potential energy, calculated from 4.0 * (1 / pow(rcut, 12) - 1/ pow(rcut, 6))
    const double rlist2;                            // square of distance threshold of neighborlist defined as 3.5^2 (reduced unit)
    Neighborlist nlist;                             // object of neighborlist
    bool need_update_nlist;                         // whether or not to update the neighborlist
    Rdf rdf;                                        // object of radial distribution function
    Property property;                              // object of mean squared displacement calculation
    double ensemble_potential;                      // potential energy of ensemble
    double ensemble_kinetic;                        // kinetic energy of ensemble
    ofstream ensemble_out;                          // output file stream of energy
    ofstream particle_out;                          // output file stream of trajectory of selected particle
    ofstream temperature_out;                       // output file stream of temperature
    ofstream msd_out;                               // output file stream of mean square displacement
};


#endif
