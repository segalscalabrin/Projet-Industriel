#ifndef INCLUDE_HPP
#define INCLUDE_HPP

#include "Neos.hpp"
#include "SimulationParameters.hpp"
#include "NeosSolution.hpp"
#include "Grid.hpp"
#include "Laplacian.hpp"

using namespace neos;

#include <cmath>
#include <iostream>
#include <cstdlib>
#include <string>
#include <vector>
#include <numeric>

struct Data {
    int dim = 2;
    
    int level = 4;
    int maxLevel = 7;

    double dt = 0.01;
    double t = 0.0;
    double tmax = 1.0;
};

#endif
