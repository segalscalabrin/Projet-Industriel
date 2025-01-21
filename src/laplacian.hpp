#ifndef LAPLACIAN_HPP
#define LAPLACIAN_HPP

#include "include.hpp"
#include "prediction-projection/Projection.hpp"

void buildLaplacianMatrix(Grid *grid, Laplacian *lap, Projection *proj, double t, PiercedVector<double> mu);

#endif
