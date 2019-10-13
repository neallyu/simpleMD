#ifndef INTEGRAL_H
#define INTEGRAL_H

#include <cmath>

using namespace std;

// gap number in integral
int n = 10000;

double lenard_jones_potential(double epsilon, double sigma, double r) {
    return 4 * epsilon * (pow(sigma, 12) / pow(r, 12) - pow(sigma, 6) / pow(r, 6));
}

// Cotes method
double integral(const double &epsilon, const double &sigma, const double &a, const double &b)
{
    double sum = 0.0;
    double gaps = (b - a) / double(n);
    double h = gaps / 2.0;
    for (int i = 0; i < n; i++)
    {
        sum += (h / 45.0) * (7.0 * lenard_jones_potential(epsilon, sigma, a + i * gaps) +
                32.0 * lenard_jones_potential(epsilon, sigma, a + i * gaps + 0.25 * gaps) +
                12.0 * lenard_jones_potential(epsilon, sigma, a + i * gaps + 0.5 * gaps) +
                32.0 * lenard_jones_potential(epsilon, sigma, a + i * gaps + 0.75 * gaps) +
                7.0 * lenard_jones_potential(epsilon, sigma, a + (i + 1) * gaps));
    }
    return sum;
}


#endif