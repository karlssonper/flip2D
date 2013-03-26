#include "multigrid.h"
#include "gaussSeidel.h"

Multigrid::Multigrid(Settings::Ptr s) : PressureSolver(s, MULTIGRID)
{
    _numFullCycles = s->numFullCycles;
    _numVCycles = s->numVCycles;
    _numPreSweeps = s->numPreSweeps;
    _numPostSweeps = s->numPostSweeps;

    //Save some space by resizing the arrays in pressureSolver;
    PressureSolver::_A.resize(0,0,1.f);

    int total = log(static_cast<float>(s->nx)) / log(2.f);
    int min = log(static_cast<float>(s->nxMin)) / log(2.f);
    _M = total - min;
    
    _fluidPhi.resize(s->nx,s->ny,_M,s->dx);
    _solidPhi.resize(s->nx,s->ny,_M,s->dx);
    _uWeights.resize(s->nx,s->ny,_M,s->dx);
    _vWeights.resize(s->nx,s->ny,_M,s->dx);
    _p.resize(s->nx,s->ny,_M,s->dx);
    _pTmp.resize(s->nx,s->ny,_M,s->dx);
    _b.resize(s->nx,s->ny,_M,s->dx);
    _r.resize(s->nx,s->ny,_M,s->dx);
    _A.resize(s->nx,s->ny,_M,s->dx);
}

void Multigrid::buildLinearSystem(Grid::Ptr grid,
                                  SolidSDF::Ptr solid,
                                  FluidSDF::Ptr fluid,
                                  float dt)
{
    // Can optimize this
    bool isSolidUpdated = true;
    if (isSolidUpdated) {
        // Copy values
        _solidPhi[_M-1].copy(solid->phi());
        _uWeights[_M-1].copy(grid->uWeights());
        _vWeights[_M-1].copy(grid->vWeights());

        //Downsample and re-create weights
        _downsampleSolidPhi();
        _createWeights();
    }


    // Copy and downsample the fluid signed distance field
    _fluidPhi[_M-1].copy(fluid->phi());
    _downsampleFluidPhi();
    _extrapolateIntoSolid();

    // Build the Laplacian linear systems
    _buildLaplacians(dt);

    // Build RHS and store in the PressureSolver b array
    // (the one in this class  will be modified)
    _buildRHS(grid->u(),
              grid->v(),
              grid->uWeights(),
              grid->vWeights(),
              fluid,
              PressureSolver::_b);
    _b[_M-1].copy(PressureSolver::_b);
}
    
void Multigrid::solveLinearSystem(FluidSDF::Ptr f, float dt)
{
    _p.reset();

    for (int m = 0; m < _numFullCycles; ++m) {
        _fullCycle();
    }

    for (int m = 0; m < _numVCycles; ++m) {
        _VCycle(m);
    }

    // Filter out prolong artifacts so theres only pressure values inside
    // the fluid
    _filter(_fluidPhi[_M-1], _p[_M-1], _pressure);
}

void Multigrid::_fullCycle()
{
    //Store the previous pressure in the PressueSolver array
    _pressure.copy(_p[_M-1]);

    _computeResidual(_fluidPhi[_M-1],_A[_M-1],_p[_M-1],_b[_M-1],_r[_M-1]);

    for (int m = _M - 2; m >= 0; --m) {
        _restrict(_r[m+1], _r[m]);
    }

    for (int m = 0; m < _M; ++m) {
        if (m) {
            _prolong(_p[m-1], _p[m]);
        } else {
            _p[m].reset();
        }

        _b[m].copy(_r[m]);
        _VCycle(m);
    }

    _p[_M-1].add(_pressure);
    _b[_M-1].copy(PressureSolver::_b);
}

void Multigrid::_VCycle(int m)
{
    _smooth(m, _numPreSweeps);
    if (m) {
        _computeResidual(_fluidPhi[m],_A[m],_p[m],_b[m],_r[m]);
        _restrict(_r[m], _b[m-1]);
        _p[m-1].reset();
        _VCycle(m-1);
        _prolongAndAdd(_p[m-1],_p[m]);
    }
    _smooth(m, _numPostSweeps);
}

void Multigrid::_smooth(int m, int iterations)
{
    for (int i = 0; i < iterations; ++i) {
        GaussSeidel::redBlackIteration(true,_fluidPhi[m],_A[m],_b[m],_p[m]);
        GaussSeidel::redBlackIteration(false,_fluidPhi[m],_A[m],_b[m],_p[m]);
    }
}

void Multigrid::_downsampleSolidPhi()
{
    for (int m = _M -2; m >= 0; --m) {
        _restrict(_solidPhi[m+1], _solidPhi[m]);
    }
}

void Multigrid::_createWeights()
{
    for (int m = _M - 2; m >= 0; --m) {
        SolidSDF::createWeights(_solidPhi[m], _uWeights[m], _vWeights[m]);
    }
}

void Multigrid::_downsampleFluidPhi()
{
    for (int m = _M -2; m >= 0; --m) {
        //For now, simply restrict. Maybe add C factor later
         _restrict(_fluidPhi[m+1], _fluidPhi[m]);
    }
}

void Multigrid::_extrapolateIntoSolid()
{
    for (int m = 0; m < _M; ++m) {
        FluidSDF::extrapolateIntoSolid(_solidPhi[m],_fluidPhi[m]);
    }
}

void Multigrid::_buildLaplacians(float dt)
{
    for (int m = 0; m < _M; ++m) {
        _buildLaplace(_uWeights[m],_vWeights[m],_fluidPhi[m],_A[m],dt);
    }
}

void Multigrid::_filter(const Array2f & phi,
                        const Array2f & source,
                        Array2f & target)
{
    for (int i = 0; i < target.nx(); ++i) {
        for (int j = 0; j < target.ny(); ++j) {
            target(i,j) = phi(i,j) < 0 ? source(i,j) : 0;
        }
    }
}

template<typename T_ARRAY>
void Multigrid::_restrict(const T_ARRAY & source, T_ARRAY & target)
{
    _bilerp<false>(source,target);
}

template<typename T_ARRAY>
void Multigrid::_prolong(const T_ARRAY & source, T_ARRAY & target)
{
    _bilerp<false>(source,target);
}

template<typename T_ARRAY>
void Multigrid::_prolongAndAdd(const T_ARRAY & source, T_ARRAY & target)
{
    _bilerp<true>(source,target);
}

template<bool T_ADD, typename T_ARRAY>
void Multigrid::_bilerp(const T_ARRAY & source, T_ARRAY & target)
{
    for (int i = 0; i < target.nx(); ++i) {
        for (int j = 0; j < target.ny(); ++j) {
            if (T_ADD) {
                target(i,j) += source.bilerp(target.pos(i,j));
            } else {
                target(i,j) = source.bilerp(target.pos(i,j));
            }
        }
    }
}

template<typename T_ARRAY>
MultigridArray<T_ARRAY>::MultigridArray(int nxMax,
                                      int nyMax,
                                      int levels,
                                      float dxMax)
{
    resize(nxMax,nyMax,levels,dxMax);
}

template<typename T_ARRAY>
void MultigridArray<T_ARRAY>::resize(int nxMax,
                                   int nyMax,
                                   int levels,
                                   float dxMax)
{
    _x.resize(levels);
    int nx = nxMax;
    int ny = nyMax;
    float dx = dxMax;
    for (int i = 0; i < levels; ++i) {
        _x[i].resize(nx, ny, dx);
        nx /= 2;
        ny /= 2;
        dx *= 2.0f;
    }
}
template<typename T_ARRAY>
void MultigridArray<T_ARRAY>::reset()
{
    for (int i = 0 ; i < _x.size(); ++i) {
        _x[i].reset();
    }
}

