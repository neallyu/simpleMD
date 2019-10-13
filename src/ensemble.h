#ifndef ENSEMBLE_H
#define ENSEMBLE_H

#include <vector>
#include <cstdlib>
#include <ctime>
#include "particle.h"
#include <iostream>

using namespace std;

class Ensemble {
public:
    Ensemble(const Particle particle) {
        ensemble.push_back(particle);
    }

    void create_ensembe(int particle_number) {
        for (int i = 0; i < particle_number; ++i) {
			srand((unsigned)time(NULL));
			cout << rand() << endl;
        }
    }



private:
    vector<Particle> ensemble;
};

#endif