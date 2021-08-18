#ifndef RADIAL_DISTRIBUTION_FUNCTION_H
#define RADIAL_DISTRIBUTION_FUNCTION_H

#include <fstream>
#include <vector>
#include <cmath>
#include <string>
#include "utils.h"

using namespace std;

class Rdf {

public:
    Rdf(int _nbins, double _BOX, string _output_path);

    void sample(vector<Particle> &ensemble);

    void normalize(int natoms);

    void output();

private:
    const int nbins;
    const double BOX;
    double bin_width;
    unsigned n;
    vector<double> g;
    ofstream rdf_output;

};

#endif