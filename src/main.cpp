#include <iostream>
#include <fstream>
#include "box.h"
#include "ensemble.h"
#include "particle.h"

using namespace std;

int main() {
    Box box1(30, 30, 30);
    Particle particle1(0.001, 0, 0, 10.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    Particle particle2(-0.001, 0, 0, 20.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    int i = 0;
    Ensemble ensemble1(3);
    while (i <= 10000) {
        if (i % 100 == 0) {
            ensemble1[2].print();
        }
        // ensemble1[1].print();
        // ensemble1[2].print();
        ensemble1.interact(ensemble1[0]);
        ensemble1.interact(ensemble1[1]);
        ensemble1.interact(ensemble1[2]);
        ensemble1[0].movement();
        ensemble1[1].movement();
        ensemble1[2].movement();
        ++i;
    }




    // ofstream fout("particle1.log");
    // while (i <= 5500000) {
    //     particle1.movement();
    //     particle2.movement();
    //     particle1.output(fout);

    //     particle1.interact(particle2);
    //     particle2.interact(particle1);
    //     if (!box1.isInBox(particle1)) {
    //         box1.rebounce(particle1);
    //     }
    //     if (!box1.isInBox(particle2)) {
    //         box1.rebounce(particle2);
    //     }
    //     ++i;
    // }


    // fout << flush;
    // fout.close();
}