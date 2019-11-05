#ifndef PROPERTY_H
#define PROPERTY_H

#include "particle.hpp"
#include <cmath>
#include <vector>

using namespace std;

class Property {

public:
    Property() { }

    void initalize(vector<Particle>& ensemble) {
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Start_status.push_back(*it);
            Velocity_0.push_back(pow(it->v_x, 2) + pow(it->v_y, 2) + pow(it->v_z, 2));
        }
    }

    double calc_mean_square_particle_displacement(vector<Particle>& ensemble) {
        double MSD(0);
        for (int i = 0; i < ensemble.size(); ++i) {
            MSD += pow((ensemble[i].pos_x - Start_status[i].pos_x), 2) + 
                pow((ensemble[i].pos_y - Start_status[i].pos_y), 2) + 
                pow((ensemble[i].pos_z - Start_status[i].pos_z), 2);
        }
        return MSD / ensemble.size();
    }

    // double calc_velocity_autocorrelation(vector<Particle>& ensemble) {
    //     double velocity_autocorr(0);
    //     for (int i = 0; i < ensemble.size(); ++i) {
    //         velocity_autocorr += (pow(ensemble[i].v_x, 2) + pow(ensemble[i].v_y, 2) + pow(ensemble[i].v_z, 2)) * Velocity_0[i];
    //     }
    //     return velocity_autocorr;
    // }

    // double calc_velocity_autocorrelation(vector<Particle>& ensemble) {
    //     return Velocity_0[1] * (pow(ensemble[1].v_x, 2) + pow(ensemble[1].v_y, 2) + pow(ensemble[1].v_z, 2));
    // }

    double calc_velocity_autocorrelation(vector<Particle> &ensemble) {
        double velocity_autocorr(0);
        for (int i = 0; i < ensemble.size(); ++i) {
            velocity_autocorr += Start_status[i].v_x * ensemble[i].v_x + Start_status[i].v_y * ensemble[i].v_y + Start_status[i].v_z * ensemble[i].v_z;
        }
        return velocity_autocorr;
    }

private:
    vector<Particle> Start_status;
    vector<double> Velocity_0;
};

#endif