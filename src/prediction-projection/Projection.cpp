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
 * @file   Projection.cpp
 * @author Antoine Gerard <antoine.gerard@inria.fr>
 * @date   Thu Oct 24 16:19:21 2019
 *
 * @brief  This file contains Projection class
 *
 * @copyright Inria
 */


#include <numeric>
#include "Projection.hpp"
#include "InterpolatorFactory.hpp"
#include "RBF.hpp"
#include "NeosAssert.hpp"
#include "Gradient.hpp"

namespace neos {

Projection::Projection(Grid* grid,
                       StencilBuilder* stencils,
                       SimulationParameters* param)
  : _grid(grid),
  _stencils(stencils),
  _hasStencil(stencils != nullptr),
  _param(param),
  _hasParam(param != nullptr),
  _psi({})
{
  lapP = LaplacianFactory::get(lapType::FINITEVOLUME, solverType::PETSC, grid);

  _mapGlobal = _grid->getCellGlobalMap();

  for (auto& cell: _grid->getCells())
  {
    _psi.emplace(cell.getId(),0.);
  }

  if (!_hasStencil) {
    _stencils = new StencilBuilder(grid);
  }
}

PiercedVector<double> Projection::FaceCenterProjection(
  const PiercedVector<darray3>& valCC,
  interpoType interpType)
{
  //FIXME: Use the buildStencil function to compute face center velocity here
  PiercedVector<double> valFC(_grid->getInterfaceCount());
  long interId;
  std::array<long, 2> owners;
  bool boundary;
  std::vector<double> valNeighs;
  int face;
  IInterpolator* interpo = InterpolatorFactory::get(interpType);
  std::vector<int> neighs, vertices;
  std::vector<NPoint> posNeighs;
  Var type;

  ASSERT(_grid->haveBoundaryConditions(), "What are velocity boundary conditions ? \n Use BoundaryCondition class to add BC and choose Grid constructor using those BC.");

  for (auto& inter: _grid->getInterfaces())
  {
    neighs.clear();
    posNeighs.clear();
    valNeighs.clear();

    interId    = inter.getId();
    owners     = inter.getOwnerNeigh();
    boundary   = inter.isBorder();
    face       = inter.getOwnerFace();

    switch(face/2)
    {
    case 0:
      type = Var::Ux;
      break;
    case 1:
      type = Var::Uy;
      break;
    default:
    // 3D need to be implemented for Var::Uz
      type = Var::Ux;
      ASSERT(false, "Unreachable statement");
      break;
    }

    valFC.emplace(interId);
    if (!boundary)
    {
      if ( _grid->getCellLevel(owners[0]) !=
           _grid->getCellLevel(owners[1]) )
      {
        neighs.push_back(owners[0]);
        posNeighs.push_back(_grid->evalCellCentroid(owners[0]));
        valNeighs.push_back(valCC[owners[0]][face/2]);

        neighs.push_back(owners[1]);
        posNeighs.push_back(_grid->evalCellCentroid(owners[1]));
        valNeighs.push_back(valCC[owners[1]][face/2]);

        // FIXME: Owners[0] is always the finest octant ?
        vertices = _grid->getFaceVertexLocalIds(face);
        for (auto &v: vertices)
        {
          auto neighBuff = _grid->findCellVertexNeighs(owners[0], v);
          for (auto &n: neighBuff)
          {
            if (std::find(neighs.begin(), neighs.end(), n)
                == neighs.end())
            {
              neighs.push_back(n);
              posNeighs.push_back(_grid->evalCellCentroid(n));
              valNeighs.push_back(valCC[n][face/2]);
            }
          }
        }

        if (interpType == interpoType::RBF)
          ((Rbf*)interpo)->setEpsilon(0.01/_grid->evalInterfaceArea(interId));

        valFC[interId] =
          interpo->computeInterpolation(_grid->evalInterfaceCentroid(interId),
                                        posNeighs,
                                        valNeighs);
      }
      else
      {
        valFC[interId] = 0.5 * valCC[owners[0]][face/2]
                         + 0.5 * valCC[owners[1]][face/2];
        continue;
      }

    }
    else
    {
      valFC[interId] = _grid->computeBoundaryValue(_grid->evalInterfaceCentroid(interId),
                                                   valCC[owners[0]][face/2],
                                                   face,
                                                   type);
    }
  }
  return valFC;
}

void Projection::ComputeProjectionStep(Solution& solNext,
                                       Solution& sol,
                                       PiercedVector<double>& PV_levelSet)
{

  lapP->zeroMatrix();
  lapP->zeroRHS();
  //Compute face center velocity
  this->computeFCNormalVelocity(solNext,
                                sol,
                                PV_levelSet);

  // Build RHS
  this->buildRHS(solNext);

  // Build Laplacian
  lapP->setStencils(_stencils);
  PiercedVector<double> kappaCC;
  PiercedVector<double> kappaFC(_grid->getInterfaceCount());

  //FIXME: Build that before
  for (auto& cell: _grid->getCells())
  {
    kappaCC.emplace(cell.getId());
    kappaCC[cell.getId()] = _grid->evalCellVolume(cell.getId());
  }
  for (auto& inter: _grid->getInterfaces())
  {
    kappaFC.emplace(inter.getId());
    kappaFC[inter.getId()] = 1.;
  }

  //FIXME: Matrix has to be build at beginning and not
  //       at each iterations
  lapP->buildFVMatrix(kappaCC,
                     kappaFC,
                     solNext.t + solNext.dt,
                     Var::P);

  // solve projection linear system: Lap(\psi) = div(u^*)
  this->solveProjection(_psi);

  // update pressure
  for (auto& cell: _grid->getCells())
  {
    solNext.Pressure[cell.getId()] = sol.Pressure[cell.getId()]
                                      + _psi[cell.getId()] * (taylorGreen_rho(PV_levelSet[cell.getId()], 1000, 1)/solNext.dt);
  }

  // Compute correction
  this->computeCorrection(solNext);
}

void Projection::computeFCNormalVelocity(Solution& solNext,
                                         Solution& sol,
                                         PiercedVector<double>& PV_levelSet)
{
  // Compute pressure gradient
  Gradient grad(_grid->getDimension(), _grid, _stencils);
  PiercedVector<std::array<double, 3> > gradPressure
    = grad.computeFVGradient(sol.Pressure,
                             sol.t,
                             Var::P);

  Var type;
  int dim(0);
  Stencil stencil;

  const PiercedVector<Stencil>& sten = _stencils->getCCToFCStencil();
  for (auto &inter: _grid->getInterfaces())
  {
    const long& interId = inter.getId();
    const std::array<long, 2>& owners = inter.getOwnerNeigh();
    const int& face = inter.getOwnerFace();
    switch(face/2)
    {
    case 0:
      type = Var::Ux;
      dim  = 0;
      break;
    case 1:
      type = Var::Uy;
      dim = 1;
      break;
    default:
    // 3D need to be implemented for Var::Uz
      type = Var::Ux;
      dim = 0;
      ASSERT(false, "Unreachable statement");
      break;
    }

    if (!inter.isBorder())
    {
      bool isUniform = ( (_grid->getCellLevel(owners[0])
                          == _grid->getCellLevel(owners[1]) )
                         && ( _grid->nbBorder(owners[0]) < 2)
                         && ( _grid->nbBorder(owners[1]) < 2) );

      stencil = sten[interId];

      if (isUniform)
        solNext.VelocityFC[interId] = 0.;
      else
      {
        solNext.VelocityFC[interId] = sol.VelocityFC[interId];
      }

      for (size_t i=0; i<stencil.neighs.size(); i++)
      {
        const long& neighId = stencil.neighs[i];
        if (!isUniform)
        {
          // U_{fc}^* =  U_{fc}^n + interp(U_{cc}^* - U_{cc}^n)
          double value = solNext.VelocityCC[neighId][dim]
                         + (solNext.dt/taylorGreen_rho(PV_levelSet[owners[0]], 1000, 1))
                         * gradPressure[neighId][dim]
                         - sol.VelocityCC[neighId][dim];

          solNext.VelocityFC[interId] +=
            stencil.weights.weights_CCToFC[i]
            * value;

        }
        else
        {
          double value = solNext.VelocityCC[neighId][dim]
                         + (solNext.dt/taylorGreen_rho(PV_levelSet[owners[0]], 1000, 1))
                         * gradPressure[neighId][dim];

          // U_{fc}^* = interp(U_{cc}^*)
          solNext.VelocityFC[interId] +=
            stencil.weights.weights_CCToFC[i]
            * value;
        }
      }
    }
    else     //Boundary interface
    {
      const NPoint& coord = _grid->evalInterfaceCentroid(interId);
      solNext.VelocityFC[interId] =
        _grid->computeBoundaryValue(coord,
                                    solNext.VelocityCC[owners[0]][dim],
                                    face,
                                    type,
                                    solNext.t+solNext.dt);
    }
  }

  const PiercedVector<double> gradFC  =
    grad.computeFCLSGradient(sol.Pressure, sol.t, Var::P);

  for (auto &inter: _grid->getInterfaces())
  {
    const long& interId = inter.getId();
    const std::array<long, 2>& owners = inter.getOwnerNeigh();

    solNext.VelocityFC[interId] -= solNext.dt/taylorGreen_rho(PV_levelSet[owners[0]], 1000, 1) * gradFC[interId];
  }
}

void Projection::buildRHS(Solution& solNext)
{
  // Compute divergence
  // FIXME: Add a method to compute divergence (Create a divergence class)
  PiercedVector<double> divergence;
  for (auto& inter:_grid->getInterfaces())
  {
    const long& interId = inter.getId();
    const auto& owners  = inter.getOwnerNeigh();
    const auto& normal = _grid->evalInterfaceNormal(interId);

    if (_grid->getCell(owners[0]).isInterior())
    {
      if (!divergence.exists(owners[0]))
        divergence.emplace(owners[0]);

      divergence[owners[0]] += solNext.VelocityFC[interId]
                               * std::accumulate(normal.begin(), normal.end(), 0)
                               * _grid->evalInterfaceArea(interId)
                               / _grid->evalCellVolume(owners[0]);

    }

    if (owners[1] != bitpit::Element::NULL_ID)
    {
      if (_grid->getCell(owners[1]).isInterior())
      {
        if (!divergence.exists(owners[1]))
          divergence.emplace(owners[1]);

        divergence[owners[1]] -= solNext.VelocityFC[interId]
                                 * std::accumulate(normal.begin(), normal.end(), 0)
                                 * _grid->evalInterfaceArea(interId)
                                 / _grid->evalCellVolume(owners[1]);

      }
    }
  }

  // Build rhs vector
  for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); cell++)
  {
    const long& cellId = cell.getId();
    int g_i = _mapGlobal.at(cellId);

    double value = divergence[cellId] * _grid->evalCellVolume(cellId);
    lapP->addRHSValue(g_i, value);
  }

}

void Projection::solveProjection(PiercedVector<double>& psi)
{
  //Solve linear system using Petsc
  lapP->toggleNeumann(true);
  lapP->solveLaplacian();
  PiercedVector<double> sol = lapP->getSolutionPV();

    for (auto cell = _grid->internalCellBegin(); cell != _grid->internalCellEnd(); cell++)
      {
        const long& id = cell.getId();

        psi[id]=sol[id];
      }
}

void Projection::computeCorrection(Solution& solNext)
{
  Gradient grad(_grid->getDimension(), _grid, _stencils);


  const auto gradFC = grad.computeFCLSGradient(_psi,
                                               solNext.t,
                                               Var::P);

  solNext.VelocityFC -= gradFC;

  // FIXME: Cell centered correction. A Face centered correction has to be
  // implemented
  const auto gradCC = grad.computeFVGradient(_psi,
                                             solNext.t,
                                             Var::P);

  solNext.VelocityCC -= gradCC;

  for (auto& cell: _grid->getCells())
  {
    const long& cellId = cell.getId();
    solNext.Ux[cellId] = solNext.VelocityCC[cellId][0];

    solNext.Uy[cellId] = solNext.VelocityCC[cellId][1];
  }

  // FIXME: Still have to had penalization here

}

// Stencil Projection::buildCCToFCStencil(const long& interId,
//                                        interpoType interpType)
// {
//
//   // Careful: stencil only built for non border interfaces
//
//   Stencil stencil;
//   std::vector<NPoint> posNeighs;
//
//   std::vector<int> vertices;
//   std::array<long, 2> owners = _grid->getInterface(interId).getOwnerNeigh();
//
//   int face = _grid->getInterface(interId).getOwnerFace();
//   if ( _grid->getCellLevel(owners[0]) !=
//        _grid->getCellLevel(owners[1]) )
//   {
//     stencil.neighs.push_back(owners[0]);
//     posNeighs.push_back(_grid->evalCellCentroid(owners[0]));
//
//
//     stencil.neighs.push_back(owners[1]);
//     posNeighs.push_back(_grid->evalCellCentroid(owners[1]));
//
//     // FIXME: Owners[0] is always the finest octant ?
//     vertices = _grid->getFaceVertexLocalIds(face);
//     for (auto &v: vertices)
//     {
//       auto neighBuff = _grid->findCellVertexNeighs(owners[0], v);
//       // FIXME: Pass vector neighs as pointer for bitpit function
//       // instead of testing with std::find
//       for (auto &n: neighBuff)
//       {
//         if (std::find(stencil.neighs.begin(), stencil.neighs.end(), n)
//             == stencil.neighs.end())
//         {
//           stencil.neighs.push_back(n);
//           posNeighs.push_back(_grid->evalCellCentroid(n));
//         }
//       }
//     }
//
//     IInterpolator* interpo = InterpolatorFactory::get(interpType);
//     if (interpType == interpoType::RBF)
//       ((Rbf*)interpo)->setEpsilon(0.01/_grid->evalInterfaceArea(interId));
//
//     stencil.weights.weights_CCToFC =
//       interpo->computeInterpolation(_grid->evalInterfaceCentroid(interId),
//                                     posNeighs);
//   }
//   else
//   {
//     stencil.neighs.push_back(owners[0]);
//     stencil.weights.weights_CCToFC.push_back(0.5);
//
//     stencil.neighs.push_back(owners[1]);
//     stencil.weights.weights_CCToFC.push_back(0.5);
//   }
//
//   return stencil;
// }
}
