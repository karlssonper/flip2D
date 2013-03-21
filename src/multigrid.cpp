#include "multigrid.h"

Multigrid::Multigrid(Settings::Ptr s) : PressureSolver(s)
{
    
}

void Multigrid::buildLinearSystem(const FaceArray2Xf & u,
                                  const FaceArray2Xf & v,
                                  const FaceArray2Xf & uWeights,
                                  const FaceArray2Xf & vWeights,
                                  FluidSDF::Ptr f,
                                  float dt)
{
    
}
    
void Multigrid::solveLinearSystem(FluidSDF::Ptr f, float dt)
{
    
}
