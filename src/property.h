#ifndef PROPERTY_H
#define PROPERTY_H

#include "particle.h"
#include "utils.h"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class Property {

public:
    Property(int _particle_number, double _time_interval, string _output_path);

    void initalize(vector<Particle>& ensemble);

    double calc_mean_square_particle_displacement(vector<Particle>& ensemble);

    // record v_x of each particle at the moment
    void sample_velocity_autocorrelation(vector<Particle> &ensemble);

    void calc_velocity_autocorrelation();

    void velocity_autocorr_output();


private:
    vector<Particle> Start_status;
    vector<vector <Particle> > Velocity;
    int v_index;
    int particle_number;
    double v_avg;
    double autocorr;
    vector<double> velocity_autocorr;
    double TIME_INTERVAL;
    ofstream velocity_autocorr_out;
};

#endif