#include "ensemble.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    cout << "[MD LOG]\tInitializing calculation..." << endl;
    // ensemble with x particles;
    Ensemble ensemble1(100, 2, 1e-5, 1e6, 20, argv[1], argv[2]);

    cout << "[MD LOG]\tStarting main interation..." << endl;
    ensemble1.iteration();
    cout << "[MD LOG]\tCalculation completed" << endl;
}
