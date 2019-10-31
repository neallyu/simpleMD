#include <fstream>
#include <string>
#include <iostream>
#include "error_handling.hpp"
#include "ensemble.hpp"
#include "utils.hpp"

using namespace std;

int main(int argc, char *argv[]) {
    string molecule;
    unsigned particle_number;
    double sigma;
    double epsilon;
    double MASS;
    double temp;
    double time_interval;
    double total_time;
    double box;

    try {
        ifstream input;
        input.open(argv[1]);
        if (!input) {
            throw ReadingFile_Open();
        }
        cout << "[MD LOG] " << get_current_time() << "\tReading the input file \"" << argv[1] << "\"..." << endl;
        int i = 0;
        string parameter;
        while (input >> parameter) {
            switch (i) {
                case 1:
                    molecule = parameter;
                    break;
                case 3:
                    particle_number = stoi(parameter);
                    break;
                case 5:
                    sigma = stod(parameter);
                    break;
                case 7:
                    epsilon = stod(parameter);
                    break;
                case 9:
                    MASS = stod(parameter);
                    break;            
                case 11:
                    temp = stod(parameter);
                    break;
                case 13:
                    time_interval = stod(parameter);
                    break;
                case 15:
                    total_time = stod(parameter);
                    break;
                case 17:
                    box = stod(parameter);
                    break;
            }
            ++i;
        }
        if (i != 18) {
            throw ReadingFile_Other();
        }
    } catch (ReadingFile_Open e) {
        cerr << "[MD ERR] "<< get_current_time() << "\tError in reading input file: " << e.what() << endl;
        return 0;
    } catch (ReadingFile_Other e) {
        cerr << "[MD ERR] "<< get_current_time() << "\tError in reading input file: " << e.what() << endl;
        return 0;
    }

    cout << "[MD LOG] " << get_current_time() << "\tParameters is successfully inputed" << endl;
    cout << "[MD LOG] " << get_current_time() << "\tmolecule: " << molecule << endl;
    cout << "[MD LOG] " << get_current_time() << "\tparticle number: " << particle_number << endl;
    cout << "[MD LOG] " << get_current_time() << "\tsigma: " << sigma << " Angstrom" << endl;
    cout << "[MD LOG] " << get_current_time() << "\tepsilon: " << epsilon << " kJ/mol"<< endl;
    cout << "[MD LOG] " << get_current_time() << "\tmass: " << MASS << " g/mol" << endl;
    cout << "[MD LOG] " << get_current_time() << "\ttemperature: " << temp << " K" << endl;
    cout << "[MD LOG] " << get_current_time() << "\ttime interval: " << time_interval << " fs" << endl;
    cout << "[MD LOG] " << get_current_time() << "\ttotal time: " << total_time << " ns" << endl;
    cout << "[MD LOG] " << get_current_time() << "\tbox size: " << box << " Angstrom" << endl;

    cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    
    // Initialize the ensemble
    Ensemble ensemble1(particle_number, sigma, epsilon, MASS, temp, time_interval, total_time, box);

    cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
    ensemble1.iteration();

    cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
}
