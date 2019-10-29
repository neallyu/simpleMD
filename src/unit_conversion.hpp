#ifndef UNIT_CONVERSION_H
#define UNIT_CONVERSION_H

#include <cmath>

class Unit_conversion {

public:
    // receive format: sigma(angstrom), epsilon(kJ/mol), mass(g/mol)
    Unit_conversion(double _sigma, double _epsilon, double _mass): 
        SIGMA(_sigma * 1e-10), EPSILON(_epsilon * 1000 / NA), MASS(_mass * 1000 / NA) { }

    double real_distance(double _distance) {
        return _distance * SIGMA;
    }

    double real_time_interval(double _time_interval) {
        return SIGMA * sqrt(MASS / EPSILON) * _time_interval;
    }

    double real_temperature(double _temperature) {
        return EPSILON * _temperature / kb;
    }

    double reduced_temerature(double _temperature) {
        return kb *_temperature / EPSILON;
    }

    double real_pressure(double _pressure) {
        return EPSILON * _pressure / pow(SIGMA, 3);
    }

    double real_density(double _density) {
        return _density / pow(SIGMA, 3);
    }

    double real_potential_energy(double _potential_energy) {
        return _potential_energy * EPSILON;
    }

    double real_kinetic_energy(double _kinetic_energy) {
        return _kinetic_energy * EPSILON;
    }

private:
    const double kb = 1.380649e-23; // bolzmann constant (J/K)
    const double NA = 6.02214076e23; // Avogadro constant

    const double SIGMA;
    const double EPSILON;
    const double MASS;
};

#endif