#include <fstream>
#include <string>
#include <iostream>
#include <sys/stat.h>
#include "ensemble.h"
#include "error_handling.h"

using namespace std;

void read_input(int argc, char *argv[]){
    unsigned particle_number;
    double init_temp;
    double set_temp;
    double time_interval;
    double equilibration_time;
    double total_time;
    double rho;

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
                    particle_number = stoi(parameter);
                    break;
                case 3:
                    init_temp = stod(parameter);
                    break;
                case 5:
                    set_temp = stod(parameter);
                    break;
                case 7:
                    time_interval = stod(parameter);
                    break;
                case 9:
                    equilibration_time = stod(parameter);
                    break;
                case 11:
                    total_time = stod(parameter);
                    break;
                case 13:
                    rho = stod(parameter);
                    break;
            }
            ++i;
        }
        if (i != 14) {
            throw ReadingFile_Other();
        }
        if (!mkdir(argv[2], S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH)) {
            throw CreatingOutputPath();
        }

    } catch (ReadingFile_Open e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
    } catch (ReadingFile_Other e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in reading input file: " << e.what() << endl;
    } catch (CreatingOutputPath e) {
        cerr << "[MD ERR] " << get_current_time() << "\tError in: " << e.what() << endl;
    }

    cout << "[MD LOG] " << get_current_time() << "\tParameters is successfully inputed" << endl;
    cout << "[MD LOG] " << get_current_time() << "\tParticle number: " << particle_number << endl;
    cout << "[MD LOG] " << get_current_time() << "\tInitial temperature: " << init_temp << endl;
    cout << "[MD LOG] " << get_current_time() << "\tEquilibration temperature: " << set_temp << endl;
    cout << "[MD LOG] " << get_current_time() << "\tEquilibration time: " << equilibration_time << endl;
    cout << "[MD LOG] " << get_current_time() << "\tTotal time: " << total_time << endl;
    cout << "[MD LOG] " << get_current_time() << "\tDensity: " << rho << endl;

    cout << "[MD LOG] " << get_current_time() << "\tInitializing calculation..." << endl;
    
    // Initialize the ensemble
    Ensemble ensemble1(particle_number, init_temp, set_temp, time_interval, equilibration_time, total_time, rho, argv[2]);

    cout << "[MD LOG] " <<  get_current_time() << "\tStarting main interation..." << endl;
    ensemble1.iteration();

    cout << "[MD LOG] " <<  get_current_time() << "\tCalculation completed" << endl;
}


int main(int argc, char *argv[]) {
    read_input(argc, argv);
}
