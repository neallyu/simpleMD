#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <random>
#include "particle.h"
#include "neighborlist.h"

using namespace std;

class Ensemble {

friend class Neighborlist;

public:
    // Initializer
    Ensemble(const unsigned, double, double, double);

    // return particle of that index
    inline Particle& operator[] (const int);
    inline void init_particle();

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
};

#endif
