#ifndef JACOBI_H_
#define JACOBI_H_

#include "pressure.h"

class Jacobi : public PressureSolver
{
  public:

    static Ptr create(Settings::Ptr s) { return new Jacobi(s); }

    virtual void solveLinearSystem(FluidSDF::Ptr fluid, float dt);

    static void iteration(const Array2f & phi,
                          const SparseLaplacianMatrix<float> & A,
                          const Array2f & b,
                          const Array2f & pFrom,
                          Array2f & p);

  protected:
    Array2f _pressureFrom;
    int _iterations;

    Jacobi(Settings::Ptr s);

    Jacobi();
    Jacobi(const Jacobi &);
    void operator&(const Jacobi&);
};

#endif
