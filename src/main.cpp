#include <iostream>
#include "box.h"

using namespace std;

int main() {
    Box box1(30, 30, 30);
    Particle particle1(0.001, 0, 0, 10.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    Particle particle2(-0.001, 0, 0, 20.0, 15.0, 15.0, 5, 10, 3.41, 0.00001);
    int i = 0;
    while (i <= 5500000) {
        particle1.movement();
        particle2.movement();
        // cout << "particle1: ";
        particle1.print_position();
        // cout << "particle2: ";
        // particle2.print_position();
        particle1.calculate_acceleration(particle2);
        particle2.calculate_acceleration(particle1);
        // particle1.interact(particle2);
        // particle2.interact(particle1);
        if (!box1.isInBox(particle1)) {
            box1.rebounce(particle1);
        }
        if (!box1.isInBox(particle2)) {
            box1.rebounce(particle2);
        }
        ++i;
    }
}