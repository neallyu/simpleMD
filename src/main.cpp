#include <iostream>
#include <fstream>
#include "box.h"
#include "ensemble.h"
#include "particle.h"

using namespace std;

int main() {
    //box size: x=30, y=30, z=30
    Box box1(nm(2), nm(2), nm(2));

    // ensemble with x particles;
    Ensemble ensemble1(5, 1e-9, 1e-3, box1);

    // define output file name
    ofstream particle_out("particle1.log");
    ofstream ensemble_out("ensemble1.log");

    // define time, which particle to show, box and output file names
    ensemble1.iteration(5e5, 1, box1, particle_out, ensemble_out);

    particle_out << flush;
    ensemble_out << flush;

    particle_out.close();
    ensemble_out.close();
}
