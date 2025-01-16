#ifndef TRANSPORT_HPP
#define TRANSPORT_HPP

#include "include.hpp"

PiercedVector<double> transportValue(Grid *grid, Transport *trpt, PiercedVector<double> &phi, PiercedVector<NPoint> vitesse, double dt);

#endif
