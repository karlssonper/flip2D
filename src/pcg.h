#ifndef PCG_H_
#define PCG_H_

#include "pressure.h"

class PCG : public PressureSolver
{
  public:
    static Ptr create(Settings::Ptr s)
    {
        return new PCG(s);
    }

    virtual void buildLinearSystem(const FaceArray2Xf & u,
                                   const FaceArray2Xf & v,
                                   const FaceArray2Xf & uWeights,
                                   const FaceArray2Xf & vWeights,
                                   FluidSDF::Ptr f,
                                   float dt);

    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt);
    
  protected:
    Array2f _z;
    Array2f _s;
    Array2f _q;
    Array2f _precon;
    float _tol;
    int _maxIterations;
    
    PCG(Settings::Ptr s);
    PCG();
    PCG(const PCG &);
    void operator=(const PCG&);

    void _applyPreconditioner(FluidSDF::Ptr f);

    void _applyLaplace(FluidSDF::Ptr f, const Array2f x, Array2f & b);
            
    void _buildIncompleteCholeskyPreconditioner(FluidSDF::Ptr f);
};

#endif
