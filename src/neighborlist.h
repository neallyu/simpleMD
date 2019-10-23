#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>
#include "particle.h"
#include "ensemble.h"

using namespace std;

class Neighborlist {

public:
    Neighborlist(vector<Particle> &_ensemble, double box, double _rlist2): 
        ensemble(_ensemble), BOX(box), rlist2(_rlist2) { }

    void update_neighbor_list(Particle &particle1, Particle&particle2) {
        #pragma omp parallel
        for (int i = 0; i < nlist.size(); i++)
        {
            nlist[i].resize(0);
        }

    // Atoms are not double counted in the neighbor list. That is, when atom j
    // is on atom i's list, the opposite is not true.
        #pragma omp parallel
        for (int i = 0; i < nlist.size() - 1; ++i) {
            for (int j = i + 1; j < nlist.size(); ++j) {
                if (calc_distance2(ensemble[i], ensemble[j]) < rlist2) {
                    nlist[i].push_back(j);
                }
            }
        }
    }

private:
    double calc_distance2(Particle& particle1, Particle& particle2) {
        double dx = particle1.pos_x - particle2.pos_x;
        double dy = particle1.pos_y - particle2.pos_y;
        double dz = particle1.pos_z - particle2.pos_z;

        // periodic boundary conditions
        dx -= BOX * round(dx / BOX);
        dy -= BOX * round(dy / BOX);
        dz -= BOX * round(dz / BOX);

        return dx * dx + dy * dy + dz * dz;
    }

    vector<vector <int> > nlist;

    vector<Particle> ensemble;

    const double BOX;

    const double rlist2;
};


#endif