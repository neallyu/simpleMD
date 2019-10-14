#ifndef BOX_H
#define BOX_H

#include <stdexcept>
#include "particle.h"
#include "ensemble.h"

class Box {
public:
    //Initializer, receive 3 arguments and create a box in the range of x(0, dim1), y(0, dim2), z(0, dim3)
    Box(double _dim1, double _dim2, double _dim3): dim1(_dim1), dim2(_dim2), dim3(_dim3) {
        if (_dim1 <=0 || _dim2 <= 0 || _dim3 <= 0) {
            throw runtime_error("Error: negative box size.");
        }
    }

    // return if the point is in the box
    bool isInBox(const Particle &particle) const {
        return particle.pos_x <= dim1 && particle.pos_y <= dim2 && particle.pos_z <= dim3 &&
            particle.pos_x >= 0 && particle.pos_y >= 0 && particle.pos_z >= 0;
    }

    // rebounce if particle hits the wall of box (particle version)
    void rebounce(Particle &particle) {
        if ((particle.pos_x >= dim1 && particle.v_x > 0) || 
            (particle.pos_x <= 0 && particle.v_x < 0) ) {
            particle.v_x = -particle.v_x;
        }
        if ((particle.pos_y >= dim2 && particle.v_y > 0) || 
            (particle.pos_y <= 0 && particle.v_y < 0) ) {
            particle.v_y = -particle.v_y;
        }
        if ((particle.pos_z >= dim3 && particle.v_z > 0) || 
            (particle.pos_z <= 0 && particle.v_z < 0) ) {
            particle.v_z = -particle.v_z;
        }
    }

    // rebounce if particle hits the wall of box (ensemble version)
    void rebounce(Ensemble &ensemble) {
        for (auto particle_ptr = ensemble.ensemble.begin(); particle_ptr != ensemble.ensemble.end(); ++particle_ptr) {
            if ((particle_ptr->pos_x >= dim1 && particle_ptr->v_x > 0) || 
                (particle_ptr->pos_x <= 0 && particle_ptr->v_x < 0) ) {
                particle_ptr->v_x = -particle_ptr->v_x;
            }
            if ((particle_ptr->pos_y >= dim2 && particle_ptr->v_y > 0) || 
                (particle_ptr->pos_y <= 0 && particle_ptr->v_y < 0) ) {
                particle_ptr->v_y = -particle_ptr->v_y;
            }
            if ((particle_ptr->pos_z >= dim3 && particle_ptr->v_z > 0) || 
                (particle_ptr->pos_z <= 0 && particle_ptr->v_z < 0) ) {
                particle_ptr->v_z = -particle_ptr->v_z;
            }
        }
    }

private:
    //box dimension
    double dim1;
    double dim2;
    double dim3;
};

#endif