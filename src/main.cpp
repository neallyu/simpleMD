#include <fstream>
#include "ensemble.hpp"
#include "particle.hpp"

using namespace std;

int main() {
    // ensemble with x particles;
    Ensemble ensemble1(30, 2, 1e-5, 20);

    // define output file name
    ofstream particle_out("particle1.log");
    ofstream ensemble_out("ensemble1.log");

    // define time, which particle to show, box and output file names
    ensemble1.iteration(5e4, 1, particle_out, ensemble_out);
 
    particle_out << flush;
    ensemble_out << flush;

    particle_out.close();
    ensemble_out.close();
}
