#ifndef PROPERTY_H
#define PROPERTY_H

#include "particle.hpp"
#include "utils.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>

using namespace std;

class Property {

public:
    Property(int _sample_point):v_index(0), autocorr(0), v_avg(0), velocity_autocorr_out("../output/velocity_autocorr.csv") { }

    void initalize(vector<Particle> ensemble) {
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Start_status.push_back(*it);
        }
    }

    double calc_mean_square_particle_displacement(vector<Particle>& ensemble) {
        double MSD(0);
        for (int i = 0; i < ensemble.size(); ++i) {
            MSD += pow((ensemble[i].pos_x - Start_status[i].pos_x), 2) + 
                pow((ensemble[i].pos_y - Start_status[i].pos_y), 2) + 
                pow((ensemble[i].pos_z - Start_status[i].pos_z), 2);
            // MSD += pow(ensemble[i].pos_x - Start_status[i].pos_x, 2);
        }
        return MSD / ensemble.size();
    }


    void sample_velocity_autocorrelation(vector<Particle> &ensemble) {
        Velocity.resize(v_index + 1);
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Velocity[v_index] += (calc_velocity(*it));
            v_avg += calc_velocity(*it);
        }
        ++v_index;
    }


    void calc_velocity_autocorrelation() {
        v_avg /= v_index;
        double v_variance(0);
        for (int i = 0; i < Velocity.size(); ++i) {
            v_variance += pow(Velocity[i] - v_avg, 2);
        }

        for (int l = 0; l <= Velocity.size(); ++l) {
            for (int i = 0; i < Velocity.size() - l; ++i) {
                autocorr += (Velocity[i] - v_avg) * (Velocity[i+l] - v_avg);
            }
            velocity_autocorr.push_back(autocorr / v_variance);
            autocorr = 0;
        }
    }


    void velocity_autocorr_output() {
        for (int i = 0; i < velocity_autocorr.size(); ++i){
            velocity_autocorr_out << i << "    " << velocity_autocorr[i] << endl;
        }
        velocity_autocorr_out.close();
    }


private:
    vector<Particle> Start_status;
    vector<double> Velocity;
    int v_index;
    double autocorr;
    double v_avg;
    vector<double> velocity_autocorr;
    ofstream velocity_autocorr_out;
};

#endif