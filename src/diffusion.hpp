#ifndef DIFFUSION_H
#define DIFFUSION_H

#include "particle.hpp"
#include <cmath>
#include <vector>

using namespace std;

class Diffusion {

public:
    Diffusion() { }

    void initalize(vector<Particle>& ensemble) {
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
        }
        return MSD;
    }

private:
    vector<Particle> Start_status;
};

#endif