#include <fstream>
#include "ensemble.hpp"
#include "particle.hpp"

using namespace std;

int main() {
    cout << "[MD LOG]\tInitializing calculation..." << endl;
    // ensemble with x particles;
    Ensemble ensemble1(100, 2, 1e-5, 20);

    cout << "[MD LOG]\tCreating output file" << endl;
    // define output file name
    ofstream particle_out("particle1.log");
    ofstream ensemble_out("ensemble1.log");

    cout << "[MD LOG]\tStarting main interation..." << endl;
    // define time, which particle to show, box and output file names
    ensemble1.iteration(5e5, 1, particle_out, ensemble_out);
 
    cout << "[MD LOG]\tCalculation completed" << endl;
    particle_out << flush;
    ensemble_out << flush;

    particle_out.close();
    ensemble_out.close();
    cout << "[MD LOG]\tOutput file saved" << endl;
}
