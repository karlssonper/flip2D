#ifndef PRESSURE_H_
#define PRESSURE_H_

#include "ptr.h"
#include "settings.h"
#include "sparse.h"
#include "sdf.h"

class PressureSolver : public SmartPtrInterface<PressureSolver>
{
  public:
    typedef SmartPtr<PressureSolver> Ptr;
    
    virtual void buildLinearSystem(const FaceArray2Xf & u,
                                   const FaceArray2Xf & v,
                                   const FaceArray2Xf & uWeights,
                                   const FaceArray2Xf & vWeights,
                                   FluidSDF::Ptr f,
                                   float dt);
    
    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt) = 0;
    
    const Array2f & pressure() const { return _pressure; }
    
  protected:
    Array2f _pressure;
    Array2f _b;
    SparseLaplacianMatrix<float> _A;
    
    PressureSolver(Settings::Ptr s);
    PressureSolver();
    PressureSolver(const PressureSolver &);
    void operator=(const PressureSolver&);

    void _buildLaplace(const FaceArray2Xf & uWeights,
                       const FaceArray2Xf & vWeights,
                       FluidSDF::Ptr f,
                       SparseLaplacianMatrix<float> & A,
                       float dt);

    void _buildRHS(const FaceArray2Xf & u,
                   const FaceArray2Xf & v,
                   const FaceArray2Xf & uWeights,
                   const FaceArray2Xf & vWeights,
                   FluidSDF::Ptr f,
                   Array2f & b);
                 
    //private:
    float _laplaceCenter(float w, float phiFluid, float phiAir);
};

#endif
