#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

class unit {

friend class Ensemble;

public:
    unit():epsilon(0), epsilon_reduced(0), sigma(0), sigma_reduced(0), r(0), r_reduced(0) {}

    void input_reduced() {

    }

    void input() {

    }

private:
    double epsilon_reduced;
    double sigma_reduced;
    double r_reduced;

    double epsilon;
    double sigma;
    double r;
    

    const double kb = 1.380649e-23; // bolzmann constant (J/K)

};

#endif