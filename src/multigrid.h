#ifndef MULTIGRID_H_
#define MULTIGRID_H_

#include "pressure.h"

class Multigrid : public PressureSolver
{
  public:
    static Ptr create(Settings::Ptr s)
    {
        return new Multigrid(s);
    }

    virtual void buildLinearSystem(const FaceArray2Xf & u,
                                   const FaceArray2Xf & v,
                                   const FaceArray2Xf & uWeights,
                                   const FaceArray2Xf & vWeights,
                                   FluidSDF::Ptr f,
                                   float dt);
    
    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt);
    
  protected:

    Multigrid(Settings::Ptr s);
    Multigrid();
    Multigrid(const Multigrid &);
    void operator=(const Multigrid &);
    
};


#endif
