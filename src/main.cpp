#include "ensemble.hpp"
#include "utils.hpp"
#include <ctime>

using namespace std;

int main(int argc, char *argv[]) {
    cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    // ensemble with x particles;
    Ensemble ensemble1(100, 2, 1e-5, 1e5, 20, argv[1], argv[2]);

    cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
    ensemble1.iteration();

    cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
}
