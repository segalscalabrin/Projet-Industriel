#include "laplacian.hpp"


void buildLaplacianMatrix(Grid *grid, Laplacian *lap, double t, PiercedVector<double> mu)
{
    std::vector<double> muInter0(0.0, grid->getInterfaces().size());
    PiercedVector<double> muInter = VtoPV(muInter0, grid);

    for(auto interface : grid->getInterfaces()) {
        long interfaceId = interface.getId();
        muInter[interfaceId] = computeInterfaceValue(grid, mu, interface);
    }

    lap->buildFVMatrix(mu, muInter, t, Ux);
}

double computeInterfaceValue(Grid *grid, PiercedVector<double> mu, bitpit::Interface interface)
{
    Polynomial2D interpolator;
    std::array<long, 3> cellNeighId = {-1, -1, -1};
    std::vector<std::array<double, 3>> cellNeighCentroid;
    std::vector<double> cellNeighValue;

    cellNeighId = getNeighId(grid, interface.getOwner());

    for(int i=0; i<3; i++) {
        cellNeighCentroid.push_back(grid->evalCellCentroid(cellNeighId[i]));
        cellNeighValue.push_back(mu[cellNeighId[i]]);
    }
    return interpolator.computeInterpolation(grid->evalInterfaceCentroid(interface.getId()), cellNeighCentroid, cellNeighValue);	
}

std::array<long, 3> getNeighId(Grid *grid, long cellId)
{
    std::array<long, 3> neighsId;
    std::vector<long> unorderedNeighs;
    unorderedNeighs = grid->findCellNeighs(cellId, 1, false);

    neighsId[2] = cellId;
    for(int i=0; i<2; i++) {
        neighsId[i] = unorderedNeighs[i];
    }

    return neighsId;
}