#include "laplacian.hpp"


void buildLaplacianMatrix(Grid *grid, Laplacian *lap, Projection *proj, double t, PiercedVector<double> mu)
{
    std::vector<double> muInter0(0.0, grid->getInterfaces().size());
    PiercedVector<double> muInter = VtoPV(muInter0, grid);

    lap->buildFVMatrix(mu, muInter, t, Ux);
}
