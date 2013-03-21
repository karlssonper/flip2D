#ifndef SDF_H_
#define SDF_H_

#include "ptr.h"
#include "array.h"
#include "particles.h"
#include "util.h"

class SolidSDF : public SmartPtrInterface<SolidSDF>
{
  public:
    typedef SmartPtr<SolidSDF> Ptr;

    static Ptr create(int nx, int ny, float dx)
    {
        return new SolidSDF(nx,ny,dx);
    }

    static float fractionInside(float phiA, float phiB);
    
    void initBoxBoundary(int width);

    float center(int i, int j) const { return _phi.center(i,j); }

    void createWeights(FaceArray2Xf & uw, FaceArray2Yf & vw);
    
  protected:
    CornerArray2f _phi;

    SolidSDF();
    SolidSDF(int nx, int ny, float dx);

    float _weight(float phiA, float phiB) const
    {
        return clamp(1.0f - fractionInside(phiA,phiB), 0.0f, 1.0f);
    }
};

class FluidSDF : public SmartPtrInterface<FluidSDF>
{
  public:
    typedef SmartPtr<FluidSDF> Ptr;

    static Ptr create(int nx, int ny, float dx)
    {
        return new FluidSDF(nx,ny,dx);
    }
    
    bool isFluid(int i, int j) const { return _phi(i,j) < 0; }

    float phi(int i, int j) const { return _phi(i,j); }
    float & phi(int i, int j) { return _phi(i,j); }

    void reconstructSurface(Particles::Ptr particles, float R, float r);

    void reinitialize(int numSwepIterations);

    void extrapolateIntoSolid(SolidSDF::Ptr solid); 
    
  protected:
    Array2f _phi;
    Array2f _sum;
    Array2<Vec2f> _pAvg;

    FluidSDF();
    FluidSDF(int nx, int ny, float dx);

    float _kernel(float s) const { return max(0.f,s*s*s*(1.0f-s*s)); }
    
    void _sweep(int i0, int i1, int j0, int j1);
    
    void _solveEikonal(float a, float b, float & phi) const;
};
    
inline void FluidSDF::_solveEikonal(float a, float b, float & phi) const
{
    const float s = sqr(a - b);
    const float p = sqrt(s) > _phi.dx() ?
            std::min(a,b) + _phi.dx() :
            0.5 * (a+b+sqrt(2*_phi.dx()*_phi.dx()-s));
    if (p < phi) {
        phi = p;
    }
}

inline float SolidSDF::fractionInside(float phiA, float phiB)  
{
    if(phiA < 0 && phiB < 0) {
        return 1;
    } else if (phiA < 0 && phiB >= 0) {
        return phiA / (phiA - phiB);
    } else if(phiA >= 0 && phiB < 0) {
        return phiB / (phiB - phiA);
    } else {
        return 0;
    }
}

#endif
