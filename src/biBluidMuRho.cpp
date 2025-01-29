#include "biBluidMuRho.hpp"

double taylorGreen_mu(double levelSet, double muF1, double muF2)
{  
  if(levelSet>=0.) {
    return muF2;
  }
  else {
    return muF1;
  }
  // return muF1 + (muF2 - muF1) / (1 + exp(-50 *(levelSet - 0.05)));
}

double taylorGreen_rho(double levelSet, double rhoF1, double rhoF2)
{
  if(levelSet>=0.) {
    return rhoF2;
  }
  else {
    return rhoF1;
  }
  // return rhoF1 + (rhoF2 - rhoF1) / (1 + exp(-50 *(levelSet - 0.05))); 
}