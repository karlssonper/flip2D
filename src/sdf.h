#ifndef SDF_H_
#define SDF_H_

#include "ptr.h"
#include "settings.h"
#include "array.h"
#include "particles.h"
#include "util.h"

template<typename T_ARRAY>
Vec2f sdfGradient(const T_ARRAY &phi, float x, float y);

class SolidSDF : public SmartPtrInterface<SolidSDF>
{
  public:
    typedef SmartPtr<SolidSDF> Ptr;

    static Ptr create(Settings::Ptr s) {return new SolidSDF(s);}

    static float fractionInside(float phiA, float phiB);

    Vec2f gradient(float x, float y) const { return sdfGradient(_phi, x, y); }

    Vec2f gradient(const Vec2f & pos) const { return gradient(pos.x,pos.y); }
    
    void initBoxBoundary(int width);

    float center(int i, int j) const { return _phi.center(i,j); }

    bool inside(const Vec2f & pos) const { return _phi.bilerp(pos) < 0; }

    void createWeights(FaceArray2Xf & uw, FaceArray2Yf & vw);

    static void createWeights(const CornerArray2f & phi,
                              FaceArray2Xf & uw,
                              FaceArray2Yf & vw);

    const CornerArray2f & phi() const { return _phi; } 
    
  protected:
    CornerArray2f _phi;
    
    SolidSDF(Settings::Ptr s);
    SolidSDF();
    SolidSDF(const SolidSDF &);
    void operator=(const SolidSDF &);

    static float _weight(float phiA, float phiB)
    {
        return clamp(1.0f - fractionInside(phiA,phiB), 0.0f, 1.0f);
    }
};

class FluidSDF : public SmartPtrInterface<FluidSDF>
{
  public:
    typedef SmartPtr<FluidSDF> Ptr;

    static Ptr create(Settings::Ptr s) {return new FluidSDF(s);}
            
    Vec2f gradient(float x, float y) const { return sdfGradient(_phi, x, y); }

    Vec2f gradient(const Vec2f & pos) const { return gradient(pos.x,pos.y); }
    
    bool isFluid(int i, int j) const { return _phi(i,j) < 0; }

    float phi(int i, int j) const { return _phi(i,j); }
    const Array2f & phi() const { return _phi; } 

    void reconstructSurface(Particles::Ptr particles, float R, float r);

    void reinitialize(int numSwepIterations);

    void extrapolateIntoSolid(SolidSDF::Ptr solid);

    static void extrapolateIntoSolid(const CornerArray2f & solid, Array2f &phi);
    
  protected:
    Array2f _phi;
    Array2f _sum;
    Array2<Vec2f> _pAvg;

    FluidSDF(Settings::Ptr s);
    FluidSDF();
    FluidSDF(const FluidSDF &);
    void operator=(const FluidSDF &);

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

template<typename T_ARRAY>
inline Vec2f sdfGradient(const T_ARRAY & phi, float x, float y)
{
    return Vec2f(phi(x+phi.dx(),y) - phi(x-phi.dx(),y),
                 phi(x,y+phi.dx()) - phi(x,y-phi.dx())) / 2.f*phi.dx();
                 
}

#endif
