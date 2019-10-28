#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

class unit {

friend class Ensemble;

public:
    // unit():epsilon(0), epsilon_reduced(0), sigma(0), sigma_reduced(0), r(0), r_reduced(0) {}

    void input_reduced() {

    }

    void input() {

    }

private:
    double epsilon_au;
    double sigma_au;
    double r_au;

    double epsilon;
    double sigma;
    double r;
    

    const double kb = 1.380649e-23; // bolzmann constant (J/K)

    const double LENGTH_UNIT = 10e-10; // m to angstrom

    const double MASS_UNIT = 1.6605391e10-24; // kg/mol to kg

    const double TIME_UNIT = 10e-12; // s to ps

};

#endif