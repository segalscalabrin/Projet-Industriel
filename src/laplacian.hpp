#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include "include.hpp"

void buildLaplacianMatrix(Grid *grid, Laplacian *lap, double t, PiercedVector<double> mu);

double computeInterfaceValue(Grid *grid, PiercedVector<double> mu, bitpit::Interface interface);

std::array<long, 3> getNeighId(Grid *grid, long cellId);

#endif
