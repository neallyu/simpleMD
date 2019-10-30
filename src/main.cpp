#include "ensemble.hpp"
#include "utils.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    // ensemble with x particles;
    // test for argon
    Ensemble ensemble1(100, 3.41, 0.996078441772, 39.948, 130, 1e-5, 1e6, 50, argv[1], argv[2]);

    cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
    ensemble1.iteration();

    cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
}
