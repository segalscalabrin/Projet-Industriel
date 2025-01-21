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

/**
 * @file   Prediction.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Mon Sep  9 13:38:26 2019
 *
 * @brief This file contains prediction class
 *
 * @copyright Inria
 */
#ifndef __PREDICTION_HPP__
#define __PREDICTION_HPP__

#include "include.hpp"

namespace neos {

class Prediction
{
public:
/**
 * @brief Constructor
 *
 * @param[in] grid: Pointer to the current grid
 * @param[in] stencils: Pointer to the StencilBuilder used to compute
 * stencils
 * @param[in] param: Pointer to the SimulationParameters.
 */
Prediction(Grid* grid,
           StencilBuilder* stencils = NULL,
           SimulationParameters* param = NULL);

/**
 * @brief destructor
 *
 */
~Prediction();

/**
 * @brief compute the prediction step in a Navier-Stokes problem
 *        See Antoine Fondaneche's Navier Stokes report.
 * @param[in]  SolPrev: Reference to the solution structure at time n-1
 * @param[in]  Sol: Reference to the solution structure at time n
 * @param[out] SolNext: Reference to the solution structure at time n+1.
 */
void ComputePredictionStep(Solution& SolPrev,
                           Solution& Sol,
                           Solution& SolNext);

/**
 * @brief Solve the prediction linear system. LHS*res = RHS.
 *
 * @param[out] res: Solution
 * @param[in]  LHS: left hand side matrix
 * @param[in]  RHS: right hand side vector
 */
void solvePrediction(PiercedVector<double> &Ux,
                     PiercedVector<double> &Uy);


private:
Grid* _grid;                   /**< Pointer to the current grid */

std::unordered_map<long, long> _mapGlobal;     /**< local to global mapping  */

SimulationParameters* _param;     /**< Pointer to simulation parameters */
bool _hasParam;                 /**<  Does the user give parameters */
Laplacian *lapU;
Laplacian *lapV;

StencilBuilder* _stencils;      /**< Pointer to the StencilBuilder */
static UserDataComm<double>* _pvComm;     /**< MPI communicator */
static int _pvComm_only_one_init_dirty;     /**< Integer to know if the
                                               communicator was initialized */
};
}

#endif
