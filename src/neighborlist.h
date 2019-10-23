#ifndef NEIGHBORLIST_H
#define NEIGHBORLIST_H

#include <vector>
#include "particle.h"
#include "ensemble.h"

using namespace std;

class Neighborlist {

friend class Ensemble;

public:
    Neighborlist(vector<Particle> &, double box, double _rlist2);

    void update_neighbor_list(vector<Particle> &);

private:
    double calc_distance2(Particle& particle1, Particle& particle2);

    vector<vector <int> > nlist;
    double BOX;
    double rlist2;
};


#endif