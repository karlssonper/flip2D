#ifndef GAUSS_SEIDEL_H_
#define GAUSS_SEIDEL_H_

#include "pressure.h"

class GaussSeidel : public PressureSolver
{
  public:

    static Ptr create(Settings::Ptr s) { return new GaussSeidel(s); }

    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt);

    static void redBlackIteration(bool red,
                                  const Array2f & phi,
                                  const SparseLaplacianMatrix<float> & A,
                                  const Array2f & b,
                                  Array2f & p);
  protected:
    int _iterations;
    
    GaussSeidel(Settings::Ptr s);
    
    GaussSeidel();
    GaussSeidel(const GaussSeidel &);
    void operator=(const GaussSeidel &);
};

#endif
