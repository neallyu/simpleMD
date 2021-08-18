#include "radial_distribution_function.h"


Rdf::Rdf(int _nbins, double _BOX, string _output_path): 
    nbins(_nbins), BOX(_BOX), n(0), 
    rdf_output(_output_path + "/rdf.csv") {
    g.resize(nbins, 0.0);
    bin_width = BOX / (2.0 * nbins);
}

void Rdf::sample(vector<Particle> &ensemble) {
    ++n;

    vector<double> g_thread(nbins, 0.0);

    for (unsigned int i = 0; i < ensemble.size() - 1; ++i) {
        for (unsigned int j = i + 1; j < ensemble.size(); ++j) {
            double d = distance(ensemble[i], ensemble[j], BOX);
            if (d < BOX / 2.0) {
                int ig = d / bin_width;
                g_thread[ig] += 2.0;
            }
        }
    }

    for (int i = 0; i < nbins; i++) {
        g[i] += g_thread[i];
    }
}


void Rdf::normalize(int natoms) {
    double norm_factor = 4.0 / 3.0 * M_PI * natoms * (natoms-1.0) * n * pow(bin_width, 3) / pow(BOX, 3);

    for (int i = 0; i < nbins; i++)
    {
        double r = (double) i;
        double binvol = pow(r+1.0, 3) - pow(r, 3);
        g[i] /= (binvol * norm_factor);
    }
}


void Rdf::output() {
    for (int i = 0; i < this->nbins; i++) {
        rdf_output << i * bin_width << "    " << g[i] << endl;
    }
    rdf_output.close();
}