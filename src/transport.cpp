#include "transport.hpp"

void transportValue(Grid *grid, Transport *trpt, PiercedVector<double> &phi, PiercedVector<NPoint> vitesse, double dt)
{
    trpt->compute(phi, vitesse, dt);
}