#ifndef MULTIGRID_H_
#define MULTIGRID_H_

#include "pressure.h"
#include <vector>

//template<typename Array>
//class MultigridArray;

template<typename T_ARRAY>
class MultigridArray
{
  public:
    MultigridArray() {}
    MultigridArray(int nxMax, int nyMax, int levels, float dxMax);
    void resize(int nxMax, int nyMax, int levels, float dxMax);
    void reset();
    const T_ARRAY & array(int level) const { return _x[level]; }
    T_ARRAY & array(int level) { return _x[level]; }
    const T_ARRAY & operator[](int level) const { return _x[level]; }
    T_ARRAY & operator[](int level) { return _x[level]; }
  protected:
    std::vector<T_ARRAY> _x;
};

class Multigrid : public PressureSolver
{
  public:
    static Ptr create(Settings::Ptr s) { return new Multigrid(s); }

    virtual void buildLinearSystem(Grid::Ptr grid,
                                   SolidSDF::Ptr solid,
                                   FluidSDF::Ptr fluid,
                                   float dt);
    
    virtual void solveLinearSystem(FluidSDF::Ptr f, float dt);
    
  protected:
    int _M;
    int _numFullCycles;
    int _numVCycles;
    int _numPreSweeps;
    int _numPostSweeps;

    MultigridArray<Array2f> _fluidPhi;
    MultigridArray<CornerArray2f> _solidPhi;
    MultigridArray<FaceArray2Xf> _uWeights;
    MultigridArray<FaceArray2Yf> _vWeights;
    MultigridArray<Array2f> _p;
    MultigridArray<Array2f> _pTmp;
    MultigridArray<Array2f> _b;
    MultigridArray<Array2f> _r;
    MultigridArray<SparseLaplacianMatrix<float> > _A;
    
    Multigrid(Settings::Ptr s);
    Multigrid();
    Multigrid(const Multigrid &);
    void operator=(const Multigrid &);

    void _fullCycle();

    void _VCycle(int m);

    void _smooth(int m, int iterations);
    
    void _downsampleSolidPhi();

    void _createWeights();

    void _downsampleFluidPhi();

    void _extrapolateIntoSolid();

    void _buildLaplacians(float dt);

    void _filter(const Array2f & phi, const Array2f & source, Array2f & target);
    
    template<typename T_ARRAY>
    void _restrict(const T_ARRAY & source, T_ARRAY & target);

    template<typename T_ARRAY>
    void _prolong(const T_ARRAY & source, T_ARRAY & target);

    template<typename T_ARRAY>
    void _prolongAndAdd(const T_ARRAY & source, T_ARRAY & target);

    template<bool T_ADD, typename T_ARRAY>
    void _bilerp(const T_ARRAY & source, T_ARRAY & target);

};




#endif
