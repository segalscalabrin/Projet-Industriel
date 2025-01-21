#include "laplacian.hpp"


void buildLaplacianMatrix(Grid *grid, Laplacian *lap, Projection *proj, double t, PiercedVector<double> mu)
{
    PiercedVector<double> muInter;

    muInter = proj->FaceCenterProjection(mu); 

    lap->buildFVMatrix(mu, muInter, t, Ux);
}
