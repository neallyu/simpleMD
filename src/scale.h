#ifndef SCALE_H
#define SCALE_H

const double kb = 1.380649e-23; // bolzmann constant (J/K)

const double AMU = 1.6605402e-27; // atom mass unit

const double NANOMETER = 10e-9;


double amu(double _amu) {
        return _amu * AMU;
    }

double nm(double _nm) {
        return _nm * NANOMETER;
    }

#endif