#ifndef GRID_H_
#define GRID_H_

#include "ptr.h"
#include "particles.h"
#include "sdf.h"

class Grid : public SmartPtrInterface<Grid>
{
  public:
    typedef SmartPtr<Grid> Ptr;

    static Ptr create(int nx, int ny, float dx) { return new Grid(nx,ny,dx); }

    void sampleVelocities(Particles::Ptr p);

    void applyGravity(const Vec2f & g, float dt);

    float CFL() const;

    void extrapolateVelocities(FluidSDF::Ptr f, int numSweepIterations);
    
  protected:
    FaceArray2Xf _u;
    FaceArray2Xf _uSum;
    FaceArray2Xf _uWeights;
    
    FaceArray2Yf _v;
    FaceArray2Yf _vSum;    
    FaceArray2Yf _vWeights;
    
    Grid(int nx, int ny, float dx);
    Grid();
    Grid(const Grid &);
    void operator=(const Grid&);

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
};

#endif
