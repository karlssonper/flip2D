#ifndef GRID_H_
#define GRID_H_

#include "ptr.h"
#include "particles.h"
#include "sdf.h"

class Grid : public SmartPtrInterface<Grid>
{
  public:
    typedef SmartPtr<Grid> Ptr;

    static Ptr create(Settings::Ptr s) { return new Grid(s); }

    void sampleVelocities(Particles::Ptr p);

    void applyGravity(const Vec2f & g, float dt);

    float CFL() const;

    void extrapolateVelocities(FluidSDF::Ptr f, int numSweepIterations);

    void pressureProjection(const Array2f &pressure, FluidSDF::Ptr f, float dt);

    void enforceBoundaryConditions(SolidSDF::Ptr s);
    
    const FaceArray2Xf & u() const { return _u;}
    const FaceArray2Xf & uWeights() const { return _uWeights;}
    const FaceArray2Yf & v() const { return _v;}
    const FaceArray2Yf & vWeigts() const { return _vWeights;}
    
  protected:
    FaceArray2Xf _u;
    FaceArray2Xf _uSum;
    FaceArray2Xf _uWeights;
    
    FaceArray2Yf _v;
    FaceArray2Yf _vSum;    
    FaceArray2Yf _vWeights;
    
    Grid(Settings::Ptr s);
    Grid();
    Grid(const Grid &);
    void operator=(const Grid&);

    template<Faces T_FACE, bool T_U>
    void _sweep(FluidSDF::Ptr f, bool upsweepX, bool upsweepY);

    template <typename T_ARRAY>
    void _accumulate(T_ARRAY & array,
                     T_ARRAY & sum,
                     T_ARRAY & marker,
                     float q,
                     int i,
                     int j,
                     float tx,
                     float ty);

    void _reset();

    bool _theta(float w, float phiA, float phiB, float & theta)
    {
        if (w > 0 && (phiA < 0 || phiB < 0)) {
            theta = 1.0f;
            if (phiA > 0 || phiB > 0) {
                theta = SolidSDF::fractionInside(phiA,phiB);
            }
            if (theta < 0.01) {
                theta = 0.01;
            }
            return true;
        } else {
            return false;
        }
    }
};

#endif
