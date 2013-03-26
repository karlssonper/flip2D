#ifndef SETTINGS_H_
#define SETTINGS_H_

#include "ptr.h"
#include "vec2.h"
#include <fstream>

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

    void write(std::ofstream & out) const
    {
        _write(out, &nx);
        _write(out, &ny);
        _write(out, &dx);
        _write(out, &initialFluidCenter);
        _write(out, &initialFluidRadius);
        _write(out, &initialVelocity);
        _write(out, &particlesPerCell);
        _write(out, &solidWidth);
        _write(out, &R);
        _write(out, &numPhiSweepIterations);
        _write(out, &gravity);
        _write(out, &numVelSweepIterations);
        _write(out, &usePCG);
        _write(out, &tolerance);
        _write(out, &maxIterations);
        _write(out, &useMultigrid);
        _write(out, &numFullCycles);
        _write(out, &numVCycles);
        _write(out, &numPreSweeps);
        _write(out, &numPostSweeps);
        _write(out, &nxMin);
        _write(out, &useGaussSeidel);
        _write(out, &numGaussSeidelIterations);
        _write(out, &useJacobi);
        _write(out, &numJacobiIterations);
    }

    void read(std::ifstream & in)
    {
        _read(in, &nx);
        _read(in, &ny);
        _read(in, &dx);
        _read(in, &initialFluidCenter);
        _read(in, &initialFluidRadius);
        _read(in, &initialVelocity);
        _read(in, &particlesPerCell);
        _read(in, &solidWidth);
        _read(in, &R);
        _read(in, &numPhiSweepIterations);
        _read(in, &gravity);
        _read(in, &numVelSweepIterations);
        _read(in, &usePCG);
        _read(in, &tolerance);
        _read(in, &maxIterations);
        _read(in, &useMultigrid);
        _read(in, &numFullCycles);
        _read(in, &numVCycles);
        _read(in, &numPreSweeps);
        _read(in, &numPostSweeps);
        _read(in, &nxMin);
        _read(in, &useGaussSeidel);
        _read(in, &numGaussSeidelIterations);
        _read(in, &useJacobi);
        _read(in, &numJacobiIterations);
    }

  protected:
    Settings() {}
    Settings(const Settings &);
    void operator=(const Settings &);
    
    template<typename T>
    void _read(std::ifstream & in, T * param)
    {
        in.read(reinterpret_cast<char *>(param), sizeof(T));
    }

    template<typename T>
    void _write(std::ofstream & out, T * param) const
    {
        out.write(reinterpret_cast<const char*>(param), sizeof(T));
    }
};

#endif
