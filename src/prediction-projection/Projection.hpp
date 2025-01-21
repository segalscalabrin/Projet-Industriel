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
 * @file   Projection.hpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Thu Aug 22 10:11:54 2019
 *
 * @brief  This file contains Projection class
 *
 * @copyright INRIA
 */

#ifndef __PROJECTION_HPP__
#define __PROJECTION_HPP__

#include "Grid.hpp"
#include "common.hpp"
#include "LaplacianFactory.hpp"
#include "NeosPiercedVector.hpp"
#include "StencilBuilder.hpp"
#include "SimulationParameters.hpp"
#include "NeosSolution.hpp"

using std::vector;

namespace neos {

class Projection
{
public:
/**
 * @brief Constructor
 *
 * @param[in] grid: Pointer to current grid
 */
Projection(Grid* grid,
           StencilBuilder* stencils = nullptr,
           SimulationParameters* param = nullptr);

/**
 * Destructor
 *
 */
~Projection(){
  ;
}

/**
 * @brief Project cell centered velocity to face center.
 *
 * @param[in] valCC: Cell centered velocity in each direction
 * @param[in] interpType: Interpolation type. Default RBF
 *
 * @return PiercedVector with normal velocity at Face center
 */
PiercedVector<double> FaceCenterProjection(
  const PiercedVector<darray3>& valCC,
  interpoType interpType = interpoType::RBF);

/**
 * @brief Compute projection step in a Navier Stokes problem
 *        See Antoine Fondaneche's team.
 * @param[in,out] solNext: reference to the solution struct at t^{n+1}
 * @param[in] sol: reference to the solution struct at t^n
 */
void ComputeProjectionStep(Solution& solNext,
                           Solution& sol);

/**
 * @brief Compute face center normal velocity for solution at time n+1
 *
 * @param[in,out] solNext: reference to the solution struct at t^{n+1}
 * @param[in] sol: reference to the solution at t^{n}
 */
void computeFCNormalVelocity(Solution& solNext,
                             Solution& sol);

/**
 * @brief Build the RHS for projection step in Navier-Stokes Problem
 *
 * @param[in] solNext:reference to the solution struct at t^{n+1}
 */
void buildRHS(Solution& solNext);

/**
 * @brief solve the linear system for the projection step
 *
 * @param[out] psi: reference to where we want to keep solution
 */
void solveProjection(PiercedVector<double>& psi);

/**
 * @brief compute the correction for face/cell-centered velocities.
 *
 * @param[in,out] solNext: reference to the solution at time t^{n+1}
 */
void computeCorrection(Solution& solNext);

/**
 * @brief Compute a stencil for Cell to Interface interpolation.
 *
 * @param interId: id of the interface we want to have stencil
 * @param interpType: type of the interpolation. Default: Radial Basis
 *                    Functions
 *
 * @return Stencil associated to the interface interId
 */
// Stencil buildCCToFCStencil(const long& interId,
//                            interpoType interpType = interpoType::RBF);

private:
Grid* _grid;                    /**< Pointer on current grid */
Laplacian *lapP;

StencilBuilder* _stencils;      /**< Pointer to the StencilBuilder */
bool _hasStencil;               /**< Does the user gives a stencil */

SimulationParameters* _param;     /**< Pointer to the SimulationParameters */
bool _hasParam;                 /**< Does the user gives built parameters */

std::unordered_map<long, long> _mapGlobal;     /**<  local to global mapping */

PiercedVector<double> _psi;     /**< solution before correction  */
};
}

#endif
