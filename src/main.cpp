#include <iostream>
#include <fstream>
#include "box.h"
#include "ensemble.h"
#include "particle.h"

using namespace std;

int main() {
    Box box1(30, 30, 30);
    // Particle particle1(0.001, 0, 0, 10.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    // Particle particle2(-0.001, 0, 0, 20.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    int i = 0;
    Ensemble ensemble1(2);
    // ofstream fout("particle1.log");
    ofstream fout2("ensemble1.log");
    while (i <= 500000) {
        ensemble1.ensemble_kinetic = 0;
        ensemble1.ensemble_potential = 0;
        // ensemble1[5].output(fout);
        for (auto particle_ptr = ensemble1.ensemble.begin(); particle_ptr != ensemble1.ensemble.end(); ++particle_ptr) {
            ensemble1.interact(*particle_ptr);
            particle_ptr->kinetic();
            ensemble1.ensemble_potential += particle_ptr->potential_value;
            ensemble1.ensemble_kinetic += particle_ptr->kinetic_value;
            particle_ptr->movement();
            box1.rebounce(ensemble1);
        }
        ensemble1.output(fout2);
        ++i;
    }
    // fout << flush;
    // fout.close();
    fout2 << flush;
    fout2.close();

}