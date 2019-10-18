#ifndef BOX_H
#define BOX_H

#include <stdexcept>

using namespace std;

class Box {

friend class Particle;

friend class Ensemble;

public:
    //Initializer, receive 3 arguments and create a box in the range of x(0, dim1), y(0, dim2), z(0, dim3)
    Box(double _dim1, double _dim2, double _dim3): dim1(_dim1), dim2(_dim2), dim3(_dim3) {
        if (_dim1 <=0 || _dim2 <= 0 || _dim3 <= 0) {
            throw runtime_error("Error: negative box size.");
        }
    }

private:
    //box dimension
    double dim1;
    double dim2;
    double dim3;
};

#endif