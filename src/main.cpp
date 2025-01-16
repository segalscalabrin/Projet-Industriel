/* -------------------------------------------------------------------------*\
 *
 *  NEOS
 *
 *  -------------------------------------------------------------------------
 *  License
 *  This file is part of Neos.
 *
 *  Neos is free software: you can redistribute it and/or modify it
 *  under the terms of the GNU Lesser General Public License v3 (LGPL)
 *  as published by the Free Software Foundation.
 *
 *  Neos is distributed in the hope that it will be useful, but WITHOUT
 *  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 *  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
 *  License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License
 *  along with Neos. If not, see <http://www.gnu.org/licenses/>.
 *
\*---------------------------------------------------------------------------*/
/* -------------------------------------------------------------------------*\
sudo docker run -it -v $(pwd):/builds/workspace neos
\*---------------------------------------------------------------------------*/


#include "include.hpp"

#include "export.hpp"
#include "transport.hpp"
#include "laplacian.hpp"


int main(int argc, char **argv) {
    Neos_Init(&argc, &argv);

    // Create struct and class pointer
    Data data;
    Grid *grid = new Grid(0.0, 0.0, 0.0, 1.0, 1.0 / std::pow(2, data.level), data.dim);
    Transport *trpt = new Transport(grid);
    Laplacian *lap = new Laplacian(grid);

    // Create pierced vector for different value
    vector<double> phi0(0.0, grid->nbCells());
    vector<double> mu0(0.0, grid->nbCells());
    vector<NPoint> vitesse0(0.0, grid->nbCells());

    PiercedVector<double> phi = VtoPV(mu0, grid);
    PiercedVector<double> mu = VtoPV(mu0, grid);
    PiercedVector<NPoint> vitesse = VtoPV(vitesse0, grid);

    while(data.t < data.tmax) {
        data.t += data.dt;
        // Compute level set 
        transportValue(grid, trpt, phi, vitesse, data.dt);

        // Compute the laplacian matrix
        buildLaplacianMatrix(grid, laplacian, data.t, mu)
    }







    delete grid;
    delete trpt;

    Neos_Finalize();

    return 0;
}


