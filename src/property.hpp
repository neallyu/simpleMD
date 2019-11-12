#ifndef PROPERTY_H
#define PROPERTY_H

#include "particle.hpp"
#include "utils.hpp"
#include <cmath>
#include <vector>
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class Property {

public:
    Property(int _sample_point, string _output_path):v_index(0), autocorr(0), v_avg(0),
        velocity_autocorr_out(_output_path + "/velocity_autocorr.csv") { }

    void initalize(vector<Particle> ensemble) {
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Start_status.push_back(*it);
        }
    }

    double calc_mean_square_particle_displacement(vector<Particle>& ensemble) {
        double MSD(0);
        for (int i = 0; i < ensemble.size(); ++i) {
            MSD += pow((ensemble[i].pos_x - Start_status[i].pos_x), 2);
        }
        return MSD / ensemble.size();
    }


    void sample_velocity_autocorrelation(vector<Particle> &ensemble) {
        Velocity.resize(v_index + 1);
        for (auto it = ensemble.begin(); it != ensemble.end(); ++it) {
            Velocity[v_index].push_back(*it);
            v_avg += it->v_x;
        }
        ++v_index;
    }


    void calc_velocity_autocorrelation() {
        cout << "\r[MD LOG] " << get_current_time() << "\t" << "Calculating velocity autocorrelation function" << endl;
        double ITERATION_PERCENTAGE(0);
        double V_VARIANCE(0);
        double particle_number = Velocity[0].size();
        v_avg /= (v_index * particle_number);
        #pragma omp parallel for schedule(dynamic)
        for (int t_frame = 0; t_frame < v_index; ++t_frame) {
            for (int particle_index = 0; particle_index < Velocity[0].size(); ++ particle_index) {
                for (int n = 0; n < v_index - t_frame; ++n) {
                    autocorr += ((Velocity[n][particle_index].v_x - v_avg) * (Velocity[n + t_frame][particle_index].v_x - v_avg));
                }
                V_VARIANCE += pow(Velocity[t_frame][particle_index].v_x - v_avg, 2);
            }
            velocity_autocorr.push_back(autocorr);
            autocorr = 0;
            ITERATION_PERCENTAGE = ((float) t_frame / (float) v_index) * 100;
            cout << "\r[MD LOG] " << get_current_time() << "\t" << ITERATION_PERCENTAGE << "\% completed " << flush;
        }
        cout << endl;
        for (auto it = velocity_autocorr.begin(); it != velocity_autocorr.end(); ++it) {
            *it /= V_VARIANCE;
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
    vector<vector <Particle> > Velocity;
    int v_index;
    double v_avg;
    double autocorr;
    vector<double> velocity_autocorr;
    ofstream velocity_autocorr_out;
};

#endif