#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "ptr.h"
#include "vec2.h"

class Settings : public SmartPtrInterface<Settings>
{
  public:
    typedef SmartPtr<Settings> Ptr;

    static Ptr create() { return new Settings(); }
    
    // Resolution and size
    int nx;
    int ny;
    float dx;

    // Particles
    Vec2f initialFluidCenter;
    float initialFluidRadius;
    Vec2f initialVelocity;
    int particlesPerCell;
    
    // SDF
    float solidWidth;
    float R;
    float r;
    int numPhiSweepIterations;

    // Grid
    Vec2f gravity;
    int numVelSweepIterations;
    
    // PCG
    bool usePCG;
    float tolerance;
    int maxIterations;

    // Multigrid
    bool useMultigrid;
    int numFullCycles;
    int numVCycles;
    int numPreSweeps;
    int numPostSweeps;
    int nxMin;

    // GaussSeidel;
    bool useGaussSeidel;
    int numGaussSeidelIterations;

    // Jacobi
    bool useJacobi;
    int numJacobiIterations;

  protected:
    Settings() {}
    Settings(const Settings &);
    void operator=(const Settings &);
};

#endif
